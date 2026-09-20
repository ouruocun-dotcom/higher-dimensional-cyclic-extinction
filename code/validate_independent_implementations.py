"""Compare the independent NumPy and Node/V8 Gillespie implementations."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
DATA_ROOT = ROOT.parent / "data"


def js(p, q):
    p, q = np.asarray(p, float), np.asarray(q, float)
    m = (p + q) / 2
    return float(sum(0.5 * np.sum(z[z > 0] * np.log(z[z > 0] / m[z > 0])) for z in (p, q)))


py = json.loads((DATA_ROOT / "independent_python_q3.json").read_text(encoding="utf-8"))
node = json.loads((DATA_ROOT / "independent_node_q3.json").read_text(encoding="utf-8"))
pyrow = py["results"][0]
n_py = py["design"]["trials"]
rng = np.random.default_rng(20261006)
out = {}

for name in ("switching", "static_naive", "static_effective_s", "static_effective_k"):
    # Reverse the Jeffreys smoothing used by the independent NumPy implementation.
    py_counts = np.rint(np.asarray(pyrow["probabilities"][name]) * (n_py + 1.5) - 0.5).astype(int)
    node_counts = np.asarray(node["designs"][name]["counts"], dtype=int)
    p_py, p_node = py_counts / py_counts.sum(), node_counts / node_counts.sum()
    observed = js(p_py, p_node)
    pooled = (py_counts + node_counts) / (py_counts.sum() + node_counts.sum())
    null = np.empty(20000)
    for b in range(len(null)):
        a = rng.multinomial(py_counts.sum(), pooled)
        c = rng.multinomial(node_counts.sum(), pooled)
        null[b] = js(a / a.sum(), c / c.sum())
    out[name] = {
        "python_counts": py_counts.tolist(),
        "node_counts": node_counts.tolist(),
        "js_between_implementations": observed,
        "monte_carlo_homogeneity_p": float((1 + np.sum(null >= observed)) / (len(null) + 1)),
    }

(DATA_ROOT / "independent_implementation_validation.json").write_text(json.dumps(out, indent=2), encoding="utf-8")
print(json.dumps(out, indent=2))
