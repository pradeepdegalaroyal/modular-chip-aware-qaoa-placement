"""Reproduce all results: python experiments/run_all.py [--quick]"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from modqaoa import (cut_edges, exact_maxcut, exact_partition, kl_partition, make_graph,  # noqa: E402
                     maxcut_diag, noisy_expectation, optimise_angles, random_partition)

P_INTRA = 0.005
KINDS = ["regular3", "community"]


def run(quick: bool) -> pd.DataFrame:
    """Sweep graph family x size x seed x error ratio x placement method."""
    seeds = range(2 if quick else 4)
    sizes = [8] if quick else [8, 10, 12]
    ratios = [1, 10] if quick else [1, 5, 20]
    traj = 40 if quick else 200
    rows = []
    for kind in KINDS:
        for n in sizes:
            for s in seeds:
                g = make_graph(kind, n, s)
                diag, cmax = maxcut_diag(g), exact_maxcut(g)
                th = optimise_angles(g, 2, diag, s, restarts=2)
                rng = np.random.default_rng(100 + s)
                parts = {"random": random_partition(n, rng), "kl": kl_partition(g, s),
                         "exact": exact_partition(g)}
                clean = noisy_expectation(g, th, diag, parts["exact"], 0.0, 0.0, 1, 0) / cmax
                for r in ratios:
                    for name, part in parts.items():
                        e = noisy_expectation(g, th, diag, part, P_INTRA, P_INTRA * r, traj, 7 + s)
                        rows.append(dict(kind=kind, n=n, seed=s, ratio=r, method=name,
                                         cut_edges=cut_edges(g, part), edges=g.number_of_edges(),
                                         ar=e / cmax))
                    rows.append(dict(kind=kind, n=n, seed=s, ratio=r, method="noiseless",
                                     cut_edges=0, edges=g.number_of_edges(), ar=clean))
                print(kind, n, s, flush=True)
    return pd.DataFrame(rows)


def figures(df: pd.DataFrame, out: Path) -> None:
    """Write the two result figures."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    out.mkdir(exist_ok=True)
    cols = {"random": "tab:red", "kl": "tab:blue", "exact": "tab:green"}
    fig, axs = plt.subplots(1, 2, figsize=(10, 4))
    for ax, kind in zip(axs, KINDS):
        n = df.n.max()
        d = df[(df.kind == kind) & (df.n == n)]
        for m, c in cols.items():
            g = d[d.method == m].groupby("ratio").ar
            ax.errorbar(g.mean().index, g.mean(), g.std(), label=m, color=c, marker="o", capsize=3)
        ax.axhline(d[d.method == "noiseless"].ar.mean(), color="k", ls="--", label="noiseless")
        ax.set_xscale("log")
        ax.set_xlabel("p_inter / p_intra")
        ax.set_ylabel("noisy approximation ratio")
        ax.set_title(f"{kind}, n={n}")
        ax.legend()
    fig.tight_layout()
    fig.savefig(out / "ratio_sweep.png", dpi=200)
    plt.close(fig)

    top = df.ratio.max()
    fig, ax = plt.subplots(figsize=(5.5, 4))
    for kind, ls in zip(KINDS, ["-", "--"]):
        d = df[(df.kind == kind) & (df.ratio == top)]
        piv = d.pivot_table(index=["n", "seed"], columns="method", values="ar").reset_index()
        gain = (piv["kl"] - piv["random"]).groupby(piv["n"])
        ax.errorbar(gain.mean().index, gain.mean(), gain.std(), ls=ls, marker="o", capsize=3, label=kind)
    ax.axhline(0, color="k", lw=0.5)
    ax.set_xlabel("qubits n")
    ax.set_ylabel(f"AR gain, KL over random (ratio {top})")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out / "scaling.png", dpi=200)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="< 1 min smoke run (writes *_quick dirs)")
    a = ap.parse_args()
    t0 = time.time()
    tag = "_quick" if a.quick else ""
    res = ROOT / f"results{tag}"
    res.mkdir(exist_ok=True)
    df = run(a.quick)
    df.to_csv(res / "raw.csv", index=False)
    s = (df.groupby(["kind", "n", "ratio", "method"])
           .agg(ar_mean=("ar", "mean"), ar_std=("ar", "std"), cut=("cut_edges", "mean"),
                edges=("edges", "mean")).reset_index())
    s.to_csv(res / "summary.csv", index=False)
    figures(df, ROOT / f"figures{tag}")
    piv = df.pivot_table(index=["kind", "n", "ratio", "seed"], columns="method", values="ar")
    piv["gain_kl"] = piv["kl"] - piv["random"]
    piv["gain_exact"] = piv["exact"] - piv["random"]
    piv["gap_exact_kl"] = piv["exact"] - piv["kl"]
    g = piv.groupby(["kind", "n", "ratio"])[["gain_kl", "gain_exact", "gap_exact_kl"]].agg(["mean", "std"])
    g.columns = ["_".join(c) for c in g.columns]
    g.reset_index().to_csv(res / "gains.csv", index=False)
    meta = dict(runtime_s=round(time.time() - t0, 1), quick=a.quick, p_intra=P_INTRA)
    (res / "meta.json").write_text(json.dumps(meta, indent=1))
    print(g.reset_index().to_string())
    print(meta)


if __name__ == "__main__":
    main()
