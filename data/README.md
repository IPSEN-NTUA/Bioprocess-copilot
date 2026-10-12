# Partner data (not in this repository)

The simulator was calibrated on real fed-batch runs from a project partner. Those files are
confidential and are **not** committed. Only the fitted kinetic parameters are part of the code
(`MODEL31_REFINED` and the tuned defaults in `simulator/simulator.py`).

The benchmark (`benchmark/`) does not need the partner data. Only the calibration and validation
scripts in `simulator/` do:

| Script | Needs |
|---|---|
| `validate_bc18.py`, `tune_bc18.py`, `compare_fit.py` | `data_for_fluxes.xlsx`, `BC18` |
| `predict_runs.py` | `data_for_fluxes.xlsx`, `BC12`, `BC15`, `BC19` |

To run them, copy those files into `data/partner/` (ignored by git) or point the scripts elsewhere:

```
export PARTNER_DATA_DIR=/path/to/partner/data
```
