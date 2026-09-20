"""Core figure for extinction-induced interaction loss and its intervention."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import binomtest


ROOT = Path(__file__).resolve().parent
REPO_ROOT = ROOT.parent if (ROOT.parent / "data").exists() else ROOT
DATA_ROOT = REPO_ROOT / "data" if (REPO_ROOT / "data").exists() else ROOT
FIGURE_ROOT = REPO_ROOT / "figures" if (REPO_ROOT / "figures").exists() else ROOT


def read(name: str):
    return json.loads((DATA_ROOT / name).read_text(encoding="utf-8"))


def groups(filename: str, q: int):
    rows = read(filename)["designs"]["switching"]["finalPairStats"]
    out = {"missing": [], "direct": [], "all": list(rows.values())}
    for key, row in rows.items():
        a, b = map(int, key.split("-"))
        direct = abs(a-b) in (1, q-1)
        out["direct" if direct else "missing"].append(row)
    ans = {}
    for label, vals in out.items():
        n = sum(v["count"] for v in vals)
        total = sum(v["sumResidenceTime"] for v in vals)
        total_sq = sum(v["sumSqResidenceTime"] for v in vals)
        mean = total/n
        var = max(0.0, (total_sq-n*mean*mean)/(n-1))
        ans[label] = {"n": n, "mean": mean, "se": np.sqrt(var/n)}
    ans["fraction_missing"] = ans["missing"]["n"] / ans["all"]["n"]
    return ans


baseline = {
    5: groups("pair_ring_q5.json", 5),
    7: groups("pair_ring_q7.json", 7),
    9: groups("leaky_ring_q9_0.json", 9),
}
q7_leak = [(0.0, "pair_ring_q7.json"), (0.02, "leaky_ring_q7_0p02.json"),
           (0.05, "leaky_ring_q7_0p05.json"), (0.1, "leaky_ring_q7_0p1.json"),
           (0.25, "leaky_ring_q7_0p25.json"), (0.5, "leaky_ring_q7_0p5.json")]
q9_leak = [(0.0, "leaky_ring_q9_0.json"), (0.1, "leaky_ring_q9_0p1.json"),
           (0.5, "leaky_ring_q9_0p5.json")]

fig, axes = plt.subplots(1, 3, figsize=(14.8, 4.25))

ax = axes[0]
qs = np.array([5,7,9])
obs = np.array([baseline[q]["fraction_missing"] for q in qs])
random_pair = (qs-3)/(qs-1)
ax.plot(qs, obs, "o-", color="#4C78A8", lw=1.8, label="observed final pair")
ax.plot(qs, random_pair, "s--", color="#9D9D9D", lw=1.5, label="uniform random pair")
ax.set_xticks(qs)
ax.set_ylim(.4,1.02)
ax.set_xlabel("initial number of species, $q$")
ax.set_ylabel("probability final pair lacks a ring edge")
ax.set_title("a  Final pairs without a nearest-neighbour edge", loc="left", fontweight="bold", fontsize=11)
ax.legend(frameon=False, fontsize=8)
ax.grid(axis="y", alpha=.2)

ax = axes[1]
for label, dx, color, marker in (("missing edge", -.12, "#B279A2", "D"), ("direct edge", .12, "#59A14F", "o")):
    key = "missing" if label.startswith("missing") else "direct"
    means=np.array([baseline[q][key]["mean"] for q in qs]);ses=np.array([baseline[q][key]["se"] for q in qs])
    ax.errorbar(qs+dx,means,yerr=1.96*ses,fmt=marker,color=color,ms=7,capsize=3,label=label)
ax.set_xticks(qs)
ax.set_xlabel("initial number of species, $q$")
ax.set_ylabel("final-pair residence time")
ax.set_title("b  Final-stage residence time by pair type", loc="left", fontweight="bold", fontsize=11)
ax.legend(frameon=False, fontsize=8)
ax.grid(axis="y", alpha=.2)

ax = axes[2]
for q, series, color, marker in ((7,q7_leak,"#54A24B","o"),(9,q9_leak,"#E45756","s")):
    xx=[];means=[];ses=[]
    for leakage, filename in series:
        row=groups(filename,q)["all"]
        xx.append(leakage);means.append(row["mean"]);ses.append(row["se"])
    means=np.asarray(means);ses=np.asarray(ses)
    ax.errorbar(xx,means/means[0],yerr=1.96*ses/means[0],color=color,marker=marker,lw=1.8,ms=6,capsize=2,label=f"q={q}")
ax.set_xlabel("strength of added non-neighbour edges")
ax.set_ylabel("normalized final residence time")
ax.set_title("c  Residence time after adding non-neighbour interactions", loc="left", fontweight="bold", fontsize=11)
ax.legend(frameon=False)
ax.grid(axis="y", alpha=.2)

fig.suptitle("Final-pair statistics in higher-dimensional cyclic systems", y=1.02)
fig.tight_layout()
FIGURE_ROOT.mkdir(parents=True, exist_ok=True)
fig.savefig(FIGURE_ROOT/"fragmentation_core_result.png",dpi=240,bbox_inches="tight")
fig.savefig(FIGURE_ROOT/"fragmentation_core_result.pdf",bbox_inches="tight")

baseline_summary={}
for q in qs:
    row=baseline[q]
    expected=(q-3)/(q-1)
    baseline_summary[str(q)]={
        **row,
        "uniform_pair_missing_probability":float(expected),
        "enrichment_binomial_p_greater":float(binomtest(row["missing"]["n"],row["all"]["n"],expected,alternative="greater").pvalue),
        "missing_over_direct_residence_ratio":float(row["missing"]["mean"]/row["direct"]["mean"]),
    }
summary={"baseline":baseline_summary,"q7_leakage":[],"q9_leakage":[]}
for q,series,key in ((7,q7_leak,"q7_leakage"),(9,q9_leak,"q9_leakage")):
    for leakage,filename in series:
        row=groups(filename,q)
        summary[key].append({"leakage":leakage,"all_mean":row["all"]["mean"],"all_se":row["all"]["se"]})
(DATA_ROOT/"fragmentation_core_summary.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
print(json.dumps(summary,indent=2))
