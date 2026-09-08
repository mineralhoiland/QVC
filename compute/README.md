# QVC compute package

Install from this directory:

```bash
pip install -e ".[dev]"
python scripts/run_verification.py
python -m pytest tests -q
```

Packages:

- `qvc` — canonical numerics (QVCCompute)
- `qvc_l1` — L1 dictionary / Lindblad / closure
- `qvc_bkt` — vortex-sector C1–C6

`qvc_l1` and `qvc_bkt` locate `qvc` in this same folder. Locked numbers:
`../locks/PARAMETERS.md`.
