"""
Minimal LLM client: chat with tool calls, and embeddings.

Two wire formats, chosen with COPILOT_LLM_API:
  ollama   native Ollama API (default; lets us set the context window, which the OpenAI-compatible
           endpoint of Ollama silently leaves at a few thousand tokens)
  openai   any OpenAI-compatible /chat/completions endpoint

Environment:
  COPILOT_LLM_API        ollama | openai                      (default ollama)
  COPILOT_LLM_BASE_URL   default http://localhost:11434       (openai: e.g. https://host/v1)
  COPILOT_LLM_API_KEY    bearer token for openai-style endpoints
  COPILOT_LLM_MODEL      default llama3.2
  COPILOT_JUDGE_MODEL    model used for grading (default: COPILOT_LLM_MODEL)
  COPILOT_EMBED_MODEL    default nomic-embed-text (always served by Ollama at COPILOT_EMBED_BASE_URL)
  COPILOT_NUM_CTX        Ollama context window, default 16384

Messages use one internal shape whatever the backend:
  {"role": "assistant", "content": str, "tool_calls": [{"id": str, "name": str, "arguments": dict}]}
  {"role": "tool", "tool_call_id": str, "name": str, "content": str}
"""
from __future__ import annotations

import json
import os
import urllib.error
import urllib.request

import numpy as np

API = os.environ.get("COPILOT_LLM_API", "ollama")
BASE_URL = os.environ.get("COPILOT_LLM_BASE_URL", "http://localhost:11434").rstrip("/")
API_KEY = os.environ.get("COPILOT_LLM_API_KEY", "")
MODEL = os.environ.get("COPILOT_LLM_MODEL", "llama3.2")
JUDGE_MODEL = os.environ.get("COPILOT_JUDGE_MODEL", MODEL)
EMBED_MODEL = os.environ.get("COPILOT_EMBED_MODEL", "nomic-embed-text")
EMBED_BASE_URL = os.environ.get("COPILOT_EMBED_BASE_URL", "http://localhost:11434").rstrip("/")
NUM_CTX = int(os.environ.get("COPILOT_NUM_CTX", "16384"))
TIMEOUT_S = int(os.environ.get("COPILOT_LLM_TIMEOUT", "600"))


class LLMError(RuntimeError):
    pass


def _post(url: str, payload: dict, key: str = "") -> dict:
    headers = {"Content-Type": "application/json"}
    if key:
        headers["Authorization"] = f"Bearer {key}"
    req = urllib.request.Request(url, data=json.dumps(payload).encode(), headers=headers)
    try:
        with urllib.request.urlopen(req, timeout=TIMEOUT_S) as r:
            return json.loads(r.read())
    except urllib.error.HTTPError as e:
        raise LLMError(f"{url}: HTTP {e.code} {e.read().decode(errors='replace')[:500]}") from e
    except (urllib.error.URLError, TimeoutError) as e:
        raise LLMError(f"{url}: {e}") from e


def _tool_specs(tools: list[dict]) -> list[dict]:
    return [{"type": "function", "function": {"name": t["name"], "description": t["description"],
                                              "parameters": t["parameters"]}} for t in tools]


def chat(messages: list[dict], tools: list[dict] | None = None, model: str | None = None,
         json_mode: bool = False) -> dict:
    """One model turn. Returns an assistant message in the internal shape."""
    model = model or MODEL
    if API == "ollama":
        wire = []
        for m in messages:
            w = {"role": m["role"], "content": m.get("content") or ""}
            if m.get("tool_calls"):
                w["tool_calls"] = [{"function": {"name": c["name"], "arguments": c["arguments"]}} for c in m["tool_calls"]]
            if m["role"] == "tool":
                w["tool_name"] = m.get("name", "")
            wire.append(w)
        payload = {"model": model, "messages": wire, "stream": False,
                   "options": {"temperature": 0, "num_ctx": NUM_CTX, "seed": 0}}
        if tools:
            payload["tools"] = _tool_specs(tools)
        if json_mode:
            payload["format"] = "json"
        msg = _post(f"{BASE_URL}/api/chat", payload)["message"]
        calls = [{"id": f"call_{i}", "name": c["function"]["name"], "arguments": c["function"].get("arguments") or {}}
                 for i, c in enumerate(msg.get("tool_calls") or [])]
    else:
        wire = []
        for m in messages:
            w = {"role": m["role"], "content": m.get("content") or ""}
            if m.get("tool_calls"):
                w["tool_calls"] = [{"id": c["id"], "type": "function",
                                    "function": {"name": c["name"], "arguments": json.dumps(c["arguments"])}}
                                   for c in m["tool_calls"]]
            if m["role"] == "tool":
                w["tool_call_id"] = m["tool_call_id"]
            wire.append(w)
        payload = {"model": model, "messages": wire, "temperature": 0}
        if tools:
            payload["tools"] = _tool_specs(tools)
        if json_mode:
            payload["response_format"] = {"type": "json_object"}
        msg = _post(f"{BASE_URL}/chat/completions", payload, API_KEY)["choices"][0]["message"]
        calls = []
        for c in msg.get("tool_calls") or []:
            try:
                args = json.loads(c["function"].get("arguments") or "{}")
            except json.JSONDecodeError:
                args = {"_unparsed": c["function"].get("arguments")}
            calls.append({"id": c["id"], "name": c["function"]["name"], "arguments": args})
    out = {"role": "assistant", "content": msg.get("content") or ""}
    if calls:
        out["tool_calls"] = calls
    return out


def chat_json(prompt: str, model: str | None = None) -> dict:
    """Single-turn call that must return a JSON object (used by the grader)."""
    text = chat([{"role": "user", "content": prompt}], model=model, json_mode=True)["content"]
    start, end = text.find("{"), text.rfind("}")
    try:
        return json.loads(text[start:end + 1])
    except (json.JSONDecodeError, ValueError) as e:
        raise LLMError(f"model did not return JSON: {text[:200]!r}") from e


def embed(texts: list[str], batch: int = 64) -> np.ndarray:
    """Unit-length embeddings from Ollama."""
    out = []
    for i in range(0, len(texts), batch):
        r = _post(f"{EMBED_BASE_URL}/api/embed", {"model": EMBED_MODEL, "input": texts[i:i + batch]})
        out.extend(r["embeddings"])
    arr = np.asarray(out, dtype=np.float32)
    return arr / np.maximum(np.linalg.norm(arr, axis=1, keepdims=True), 1e-9)
