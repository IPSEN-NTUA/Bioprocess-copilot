"""
Retrieval index over the raw corpus files (the RAG arm).

Every file is turned into text chunks that keep the file's own wording, headers and number formats:
nothing is normalised here, by design. Tabular files are cut into row windows with the header repeated.
Retrieval is dense (Ollama embeddings), lexical (BM25) or both interleaved (default).

Usage:  python -m copilot.rag            builds and caches the index
"""
from __future__ import annotations

import hashlib
import json
import math
import re
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import openpyxl

from copilot import llm
from copilot.build_kg import CORPUS, ROOT

CACHE = ROOT / "copilot" / ".cache"
ROWS_PER_CHUNK = 12            # 2 h of a 10-min log
MAX_CHARS = 1800


@dataclass
class Chunk:
    source: str                # corpus-relative path, plus sheet / row range
    text: str


# ---------------------------------------------------------------------------
# Chunking
# ---------------------------------------------------------------------------
def _windows(header: list[str], rows: list[str], source: str, size: int) -> list[Chunk]:
    head = "\n".join(header)
    return [Chunk(f"{source} rows {i + 1}-{min(i + size, len(rows))}", f"[{source}]\n{head}\n" + "\n".join(rows[i:i + size]))
            for i in range(0, len(rows), size)]


def _grouped(header: list[str], rows: list[str], keys: list[str], source: str) -> list[Chunk]:
    """Keep the rows of one key (e.g. one sample) together."""
    out, start = [], 0
    for i in range(1, len(rows) + 1):
        if i == len(rows) or keys[i] != keys[start]:
            out.append(Chunk(f"{source} {keys[start]}", f"[{source}]\n" + "\n".join(header) + "\n" + "\n".join(rows[start:i])))
            start = i
    return out


def _chunk_markdown(text: str, source: str) -> list[Chunk]:
    if len(text) <= MAX_CHARS:
        return [Chunk(source, f"[{source}]\n{text}")]
    title = text.split("\n", 1)[0]
    parts = re.split(r"\n(?=## )", text)
    return [Chunk(f"{source} #{i}", f"[{source}]\n{title}\n{p}" if i else f"[{source}]\n{p}") for i, p in enumerate(parts)]


def _chunk_delimited(path: Path, source: str) -> list[Chunk]:
    lines = path.read_text(encoding="utf-8").replace("\r\n", "\n").splitlines()
    n_head = next(i for i, ln in enumerate(lines) if re.match(r"^\d", ln) or (i > 1 and "," in ln and not ln.startswith("#")
                                                                              and not ln.lower().startswith(("timestamp", "inj", "sample_id", "date"))))
    header, rows = lines[:n_head], lines[n_head:]
    if "hplc/raw" in source:
        return _grouped(header, rows, [r.split(",")[2] for r in rows], source)
    size = 25 if "lims" in source or "maintenance" in source else ROWS_PER_CHUNK
    return _windows(header, rows, source, size)


def _cell(v) -> str:
    return "" if v is None else str(v)


def _chunk_workbook(path: Path, source: str) -> list[Chunk]:
    out = []
    wb = openpyxl.load_workbook(path, data_only=True, read_only=True)
    for ws in wb.worksheets:
        if ws.title == "About":
            continue
        rows = [" | ".join(_cell(v) for v in r).rstrip(" |") for r in ws.iter_rows(values_only=True)]
        src = f"{source} sheet '{ws.title}'"
        if ws.title == "Results":                                   # HPLC results: one chunk per sample
            body = [r for r in rows[1:] if r]
            out += _grouped(rows[:1], body, [r.split(" | ")[0] for r in body], src)
        elif ws.title == "GA4":
            out += _windows(rows[:5], rows[5:], src, ROWS_PER_CHUNK)
        elif len(rows) > 30:
            out += _windows(rows[:1], rows[1:], src, 20)
        else:
            out.append(Chunk(src, f"[{src}]\n" + "\n".join(rows)))
    return out


def chunk_corpus(corpus: Path = CORPUS) -> list[Chunk]:
    if "answer_key" in corpus.parts:
        raise ValueError("The answer key must never be indexed.")
    chunks = []
    for path in sorted(p for p in corpus.rglob("*") if p.is_file()):
        source = path.relative_to(corpus).as_posix()
        if path.suffix == ".md":
            chunks += _chunk_markdown(path.read_text(encoding="utf-8"), source)
        elif path.suffix == ".xlsx":
            chunks += _chunk_workbook(path, source)
        else:
            chunks += _chunk_delimited(path, source)
    return chunks


# ---------------------------------------------------------------------------
# Index
# ---------------------------------------------------------------------------
def _tokens(text: str) -> list[str]:
    return re.findall(r"[a-z0-9]+(?:[-.][a-z0-9]+)*", text.lower())


class Index:
    def __init__(self, chunks: list[Chunk], embeddings: np.ndarray | None):
        self.chunks, self.emb = chunks, embeddings
        docs = [_tokens(c.text) for c in chunks]
        self.tf = [Counter(d) for d in docs]
        self.len = np.array([len(d) for d in docs], dtype=np.float32)
        df = Counter(t for d in self.tf for t in d)
        n = len(chunks)
        self.idf = {t: math.log(1 + (n - k + 0.5) / (k + 0.5)) for t, k in df.items()}
        self.postings: dict[str, list[int]] = {}
        for i, d in enumerate(self.tf):
            for t in d:
                self.postings.setdefault(t, []).append(i)

    def _bm25(self, query: str, k1=1.5, b=0.75) -> np.ndarray:
        scores = np.zeros(len(self.chunks), dtype=np.float32)
        avg = float(self.len.mean())
        for t in set(_tokens(query)):
            for i in self.postings.get(t, []):
                f = self.tf[i][t]
                scores[i] += self.idf[t] * f * (k1 + 1) / (f + k1 * (1 - b + b * self.len[i] / avg))
        return scores

    def search(self, query: str, k: int = 8, mode: str = "hybrid", max_per_file: int = 3) -> list[Chunk]:
        """Top-k chunks. At most max_per_file come from one file, so that a long log, whose every chunk repeats
        the same header, cannot fill the whole result."""
        ranks = []
        if mode in ("dense", "hybrid"):
            if self.emb is None:
                raise RuntimeError("index was built without embeddings")
            q = llm.embed([f"search_query: {query}"])[0]
            ranks.append(np.argsort(-(self.emb @ q)))
        if mode in ("lexical", "hybrid"):
            ranks.append(np.argsort(-self._bm25(query), kind="stable"))
        # Interleave the rankings (dense 1st, lexical 1st, dense 2nd, ...). Reciprocal-rank fusion was tried first
        # and dropped: it buries a chunk that only one ranking likes, and here the lexical ranking is often the only
        # one that can match an identifier such as FB006-S08.
        merged = [i for pair in zip(*ranks) for i in pair]
        out, per_file = [], Counter()
        seen = set()
        for i in merged:
            if i in seen:
                continue
            seen.add(i)
            name = re.split(r" (?:sheet|rows|#)", self.chunks[i].source)[0].split(" ")[0]
            if per_file[name] < max_per_file:
                per_file[name] += 1
                out.append(self.chunks[i])
                if len(out) == k:
                    break
        return out


def load_index(corpus: Path = CORPUS, embed: bool = True) -> Index:
    """Build the index, or load it from the cache when the corpus and embedding model are unchanged."""
    chunks = chunk_corpus(corpus)
    if not embed:
        return Index(chunks, None)
    digest = hashlib.sha256((llm.EMBED_MODEL + "\x00".join(c.text for c in chunks)).encode()).hexdigest()[:16]
    cache = CACHE / f"rag_{digest}.npy"
    if cache.exists():
        return Index(chunks, np.load(cache))
    emb = llm.embed([f"search_document: {c.text}" for c in chunks])
    CACHE.mkdir(parents=True, exist_ok=True)
    np.save(cache, emb)
    (CACHE / f"rag_{digest}.json").write_text(json.dumps([asdict(c) for c in chunks]))
    return Index(chunks, emb)


if __name__ == "__main__":
    idx = load_index()
    by_folder = Counter(c.source.split("/")[0] for c in idx.chunks)
    print(f"{len(idx.chunks)} chunks, embeddings {idx.emb.shape}")
    for folder, n in sorted(by_folder.items()):
        print(f"  {folder:<14} {n}")
