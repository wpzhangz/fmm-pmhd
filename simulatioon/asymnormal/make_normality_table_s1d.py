#!/usr/bin/env python3
"""Generate the Theorem 3(i) coverage table for the well-conditioned
Scenario 1(d) from a refit-only run, validating the plug-in asymptotic
variance J_S^{-1} directly (no bootstrap).

Usage:
    python3 make_normality_table_s1d.py <run_dir> [--bwc=2.0] [--out=normality_table_s1d.tex]

<run_dir> must contain coverage_onestep_vs_debiased.csv and run_diagnostics.csv
produced by scenario1c_debiased_normality_implementable.R --refit-only --scenario=s1d.
"""
import csv, os, sys

TEX = {"pi": r"\pi", "mu": r"\mu", "sigma": r"\sigma",
       "alpha": r"\alpha", "tau": r"\tau"}
ORDER = ["pi1", "pi2", "mu1", "mu2", "mu3", "sigma1", "sigma2", "sigma3",
         "alpha1", "alpha2", "alpha3", "tau1", "tau2", "tau3"]


def tex_name(p):
    base = p.rstrip("0123456789")
    return "$%s_%s$" % (TEX[base], p[len(base):])


def main():
    args = sys.argv[1:]
    if not args:
        sys.exit(__doc__)
    run_dir = args[0]
    bwc = 2.0
    out = "normality_table_s1d.tex"
    for a in args[1:]:
        if a.startswith("--bwc="):
            bwc = float(a.split("=", 1)[1])
        elif a.startswith("--out="):
            out = a.split("=", 1)[1]

    cov_csv = os.path.join(run_dir, "normality_coverage.csv")
    if not os.path.exists(cov_csv):  # backward compatibility with older runs
        cov_csv = os.path.join(run_dir, "coverage_onestep_vs_debiased.csv")
    with open(cov_csv) as f:
        rows = list(csv.DictReader(f))
    with open(os.path.join(run_dir, "run_diagnostics.csv")) as f:
        diag = {int(float(r["n"])): r for r in csv.DictReader(f)}

    ns = sorted({int(float(r["n"])) for r in rows})
    d = {(int(float(r["n"])), r["parameter"]): r for r in rows}
    nb = len(ns)

    rate = ", ".join("%.1f\\%%" % (100 * float(diag[n]["success_rate"])) for n in ns)
    nlist = ", ".join(str(n) for n in ns)
    bw = "%g" % bwc

    L = [r"\begin{table}[t]", r"\centering", r"\small"]
    L.append(r"""\caption{Gaussian approximation of Theorem~\ref{thm:normality}\textup{(i)} under the
well-conditioned Scenario~1(d), which keeps the locations and weights of Scenario~1(c) but uses
narrower components and shapes away from the cusp, so that
$n\lambda_{\min}(\mathcal J_{\mathcal S})\approx3.8,7.6,15.3$ at $n=%s$ ($B=200$ replications, inference
bandwidth $h_n=%s\,n^{-1/3}$). Here $\bar z=\sqrt n\,(\widehat\theta-\theta^*)/\mathrm{SD}$,
$\widehat{\mathrm{SD}}$ is the empirical standard deviation, $\mathrm{SD}=\{\mathrm{diag}
(\mathcal J_{\mathcal S}^{-1})\}^{1/2}/\sqrt n$ the plug-in value, and Cov the empirical $95\%%$
coverage of $\widehat\theta\pm z_{0.975}\,\mathrm{SD}$. Empirical and theoretical standard deviations
agree across all coordinates and coverage is close to nominal and stable in $n$; the residual
$O(h_n^{2\underline a})$ smoothing bias leaves a mild positive shift in the scale and outer-shape
parameters that decreases with $n$.}""" % (nlist, bw))
    L.append(r"\label{tab:normality-coverage}")
    L.append(r"\begin{tabular}{l" + "cccc" * nb + "}")
    L.append(r"\toprule")
    L.append(" & " + " & ".join(r"\multicolumn{4}{c}{$n=%d$}" % n for n in ns) + r" \\")
    L.append(" ".join(r"\cmidrule(lr){%d-%d}" % (2 + 4 * i, 5 + 4 * i) for i in range(nb)))
    hdr = r"$\bar z$ & $\widehat{\mathrm{SD}}$ & $\mathrm{SD}$ & Cov"
    L.append("Parameter & " + " & ".join([hdr] * nb) + r" \\")
    L.append(r"\midrule")

    for p in ORDER:
        cells = []
        for n in ns:
            r = d[(n, p)]
            cells.append("%+.2f & %.4f & %.4f & %.3f" % (
                float(r["mean_init"]), float(r["empirical_sd"]),
                float(r["theoretical_sd"]), float(r["cov95_init"])))
        L.append(tex_name(p) + " & " + " & ".join(cells) + r" \\")

    L.append(r"\midrule")
    mc = []
    for n in ns:
        m = sum(float(d[(n, p)]["cov95_init"]) for p in ORDER) / len(ORDER)
        mc.append(r"\multicolumn{3}{c}{} & %.3f" % m)
    L.append(r"\textit{mean coverage} & " + " & ".join(mc) + r" \\")
    L.append(r"\bottomrule")
    L.append(r"\end{tabular}")
    L.append(r"\end{table}")

    path = out if os.path.dirname(out) else os.path.join(run_dir, out)
    with open(path, "w") as f:
        f.write("\n".join(L) + "\n")
    print("written:", path)


if __name__ == "__main__":
    main()
