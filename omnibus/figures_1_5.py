"""Replot Tasks 1--5 from archived numerical evidence, without changing it.

The only new scientific table is a transparent Pareto calculation from Task 3
QSSA efficiency and the archived forward GSH half-life. No outcome is imputed.
This module deliberately needs neither seaborn nor an RDKit recomputation.
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.lines import Line2D
import numpy as np

COLORS = ["#226DA8", "#168A82", "#DC863A", "#8261A8", "#B8596B", "#577648", "#645947", "#C175B5"]
INK = "#18364D"
STYLE = {
    "font.family": "DejaVu Sans", "font.size": 10,
    "axes.titlesize": 11, "axes.labelsize": 10,
    "axes.titleweight": "bold", "axes.titlecolor": INK,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.alpha": .20, "grid.linewidth": .65,
    "axes.axisbelow": True, "legend.frameon": False, "legend.fontsize": 8,
    "xtick.labelsize": 9, "ytick.labelsize": 9,
    "figure.facecolor": "white", "savefig.facecolor": "white",
    "lines.linewidth": 1.9,
}


def rows(repo: Path, relative: str) -> list[dict]:
    with (repo / relative).open(encoding="utf-8-sig", newline="") as handle:
        result = list(csv.DictReader(handle))
    if not result:
        raise ValueError(f"Empty source: {relative}")
    return result


def vals(data: list[dict], key: str) -> np.ndarray:
    result = np.array([float(row[key]) for row in data], dtype=float)
    if not np.isfinite(result).all():
        raise ValueError(f"Missing/nonfinite plotted values: {key}")
    return result


def select(data: list[dict], **terms) -> list[dict]:
    selected = [row for row in data if all(
        float(row[key]) == value if isinstance(value, (float, int)) else row[key] == value
        for key, value in terms.items())]
    if not selected:
        raise ValueError(f"Requested selection is absent: {terms}")
    return selected


def grid(data: list[dict], xkey: str, ykey: str, zkey: str):
    """Build a complete rectangular grid; absent cells are errors, never zeros."""
    xx = np.unique(vals(data, xkey))
    yy = np.unique(vals(data, ykey))
    zz = np.full((len(yy), len(xx)), np.nan)
    xlookup, ylookup = dict(zip(xx, range(len(xx)))), dict(zip(yy, range(len(yy))))
    for row in data:
        ix, iy = xlookup[float(row[xkey])], ylookup[float(row[ykey])]
        if np.isfinite(zz[iy, ix]):
            raise ValueError(f"Duplicate grid coordinate: {xkey}, {ykey}")
        zz[iy, ix] = float(row[zkey])
    if not np.isfinite(zz).all():
        raise ValueError(f"Incomplete grid: {zkey}")
    return xx, yy, zz


def canvas(task: int, title: str, footer: str):
    fig = plt.figure(figsize=(15.6, 5.7))
    gs = fig.add_gridspec(1, 3, left=.055, right=.975, bottom=.24,
                          top=.77, wspace=.35)
    fig.text(.045, .95, f"AI4PHARM  /  TASK {task:02d}", color=COLORS[0],
             fontsize=10, weight="bold")
    fig.text(.045, .885, title, color=INK, fontsize=19, weight="bold")
    fig.text(.045, .055, footer, color="#465B6A", fontsize=9, linespacing=1.55)
    return fig, gs


def save(fig, out: Path, name: str):
    fig.savefig(out / name, dpi=300, metadata={"Software": "AI4Pharm omnibus numerical evidence atlas"})
    plt.close(fig)


def task1(repo: Path, out: Path) -> dict:
    source = "projects/task01_lead_developability/results_task1/developability_results.csv"
    data = rows(repo, source)
    groups = ["Oral Drugs", "Toxic Dropouts", "bRo5 Modalities"]
    fig, gs = canvas(1, "Developability across 30 reference compounds",
        "Evidence: computed molecular descriptors + assumed ionization + unvalidated ADMET surrogates; cohort labels are descriptive.\n"
        "Radar: fixed 0--1 desirability anchors, higher is favored; the hERG axis covers one liability only. bRo5 membership does not establish oral success.")
    ax = fig.add_subplot(gs[0])
    grouped = [vals(select(data, archetype=g), "cns_mpo_approx") for g in groups]
    violin = ax.violinplot(grouped, positions=[1, 2, 3], showextrema=False, widths=.78)
    for body, color in zip(violin["bodies"], COLORS):
        body.set(facecolor=color, edgecolor=color, alpha=.24)
    for i, (points, color) in enumerate(zip(grouped, COLORS), 1):
        # Fixed positions preserve reproducibility without implying independent repeats.
        ax.scatter(i + np.linspace(-.10, .10, len(points)), np.sort(points), color=color, s=20, zorder=3)
        ax.plot([i-.17, i+.17], [np.median(points)]*2, color=INK, lw=2.5)
    ax.set(title="A  |  Approximate CNS-MPO", ylabel="Continuous desirability sum (0--6)", ylim=(-.1, 6.2))
    ax.set_xticks([1, 2, 3], ["Oral\nreferences", "Dropout\nreferences", "bRo5\nmodalities"])
    ax.text(.03, .05, "n = 10 per group; bar = median", transform=ax.transAxes, fontsize=8, color=INK)

    ax = fig.add_subplot(gs[1], projection="polar")
    angles = np.linspace(0, 2*np.pi, 6, endpoint=False)
    labels = ["ESOL", "Permeability*", "Low hERG\nliability*", "Polarity", r"Fsp$^3$", "CNS-MPO"]
    for name, color in zip(["Diazepam", "Terfenadine", "Venetoclax"], COLORS):
        row = select(data, name=name)[0]
        v = [float(row[k]) for k in ["oral_d_solubility", "oral_d_permeability", "oral_d_herg_proxy", "oral_d_polarity", "oral_d_fsp3"]]
        v.append(float(row["cns_mpo_approx"])/6)
        ax.plot(np.r_[angles, angles[0]], np.r_[v, v[0]], color=color, label=name)
        ax.fill(np.r_[angles, angles[0]], np.r_[v, v[0]], color=color, alpha=.055)
    ax.set_theta_offset(np.pi/2)
    ax.set_theta_direction(-1)
    ax.set_xticks(angles, labels, fontsize=9)
    ax.tick_params(axis="x", pad=8)
    ax.set_ylim(0, 1)
    ax.set_yticks([.25, .5, .75, 1], ["", ".5", "", "1"], fontsize=7)
    ax.set_title("B  |  Six-axis proxy profiles", pad=30)
    ax.legend(loc="upper center", bbox_to_anchor=(.5, -.12), ncol=3, columnspacing=.8, fontsize=8)

    ax = fig.add_subplot(gs[2])
    ax.axhline(500, ls="--", color="#97A5AD", lw=1)
    ax.axvline(5, ls="--", color="#97A5AD", lw=1)
    for group, color in zip(groups, COLORS):
        sub = select(data, archetype=group)
        ax.scatter(vals(sub, "clogp_rdkit"), vals(sub, "mw_da"),
                   s=24+90*vals(sub, "fsp3"), alpha=.8, color=color, edgecolor="white", linewidth=.6)
    for name, offset in [("Venetoclax", (-61, -15)), ("Cyclosporine A", (-80, -14)), ("ARV-110", (6, 7))]:
        row = select(data, name=name)[0]
        ax.annotate(name, (float(row["clogp_rdkit"]), float(row["mw_da"])),
                    xytext=offset, textcoords="offset points", fontsize=8)
    ax.set(title="C  |  Chemical space and bRo5", xlabel="RDKit cLogP", ylabel="Molecular weight (Da)")
    ax.set_ylim(0, max(vals(data, "mw_da"))*1.17)
    ax.text(.03, .96, r"Point area increases with Fsp$^3$", transform=ax.transAxes, va="top", fontsize=8)
    name = "fig_task1_developability_mpo.png"
    save(fig, out, name)
    return dict(task=1, file=name, sources=[source], panels=["CNS-MPO violin and all 30 compounds", "Six fixed-anchor desirability axes for three named references", "Molecular weight versus RDKit cLogP, area by Fsp3"],
                caveats=["Approximate ionization and unvalidated ADMET proxies", "Descriptive selected cohorts, not a clinical classifier", "No chameleonic adaptation inferred from these descriptors"])


def task2(repo: Path, out: Path) -> dict:
    base = "projects/task02_tpd/data_task2/"
    source = [base+n for n in ["equilibrium_species.csv", "affinity_cooperativity_heatmap.csv", "linker_conformers.csv", "equilibrium_metrics.csv"]]
    eq, heat, link, metrics = [rows(repo, path) for path in source]
    fig, gs = canvas(2, "Ternary binding, cooperativity and linker geometry",
        "Evidence: finite-total mass-action equilibrium and unweighted ETKDG/MMFF conformer ensembles; linker geometry is not a binding free energy.\n"
        "Hook curves are conditional on the modeled system. Their peak is an in-vitro concentration, not a clinical dose; the supplied approximate optimum is not substituted.")
    ax = fig.add_subplot(gs[0])
    for alpha, color in zip([.01, 1., 10., 100.], COLORS):
        sub = select(eq, alpha=alpha)
        ax.plot(vals(sub, "P_total_M"), vals(sub, "EPT_nM"), color=color, label=fr"$\alpha$ = {alpha:g}")
    peak = float(metrics[0]["exact_peak_total_p_nm"])
    ax.axvline(peak*1e-9, color=INK, ls=":", lw=1)
    ax.set(xscale="log", yscale="log", title="A  |  Nine-order hook curves", xlabel="Total PROTAC concentration (M)", ylabel="Ternary complex [EPT] (nM)", xlim=(1e-12, 1e-3))
    ax.legend(loc="lower center", ncol=2)
    ax.text(.02, .97, f"Exact peak: {peak:.1f} nM", transform=ax.transAxes, va="top", fontsize=8)
    ax = fig.add_subplot(gs[1])
    xx, yy, zz = grid(select(heat, alpha=10.), "Kd_E_nM", "Kd_T_nM", "peak_EPT_nM")
    im = ax.pcolormesh(xx, yy, zz, shading="auto", cmap="viridis", vmin=0)
    ax.set(xscale="log", yscale="log", title=r"B  |  Affinity sensitivity ($\alpha$ = 10)", xlabel=r"E3 binary $K_D$ (nM)", ylabel=r"Target binary $K_D$ (nM)")
    cb = fig.colorbar(im, ax=ax, fraction=.046, pad=.025)
    cb.set_label("Peak [EPT] (nM)", fontsize=9)
    ax = fig.add_subplot(gs[2])
    bins = np.linspace(min(vals(link, "distance_A"))-.1, max(vals(link, "distance_A"))+.1, 27)
    for linker, color in zip(["Flexible PEG", "Rigid alkynyl"], COLORS):
        sub = select(link, linker=linker)
        ax.hist(vals(sub, "distance_A"), bins=bins, color=color, alpha=.63,
                label=f"{linker} (n={len(sub)})", edgecolor="white", linewidth=.6)
    ax.set(title="C  |  Exit-vector distance ensemble", xlabel=r"End-to-end distance $r_{E-T}$ ($\AA$)", ylabel="Conformer count")
    ax.legend(loc="upper left")
    ax.text(.03, .60, "Unweighted counts\nRigid ensemble overlaps\nin one narrow bin", transform=ax.transAxes, fontsize=8, color=INK)
    name = "fig_task2_tpd_hook_effect.png"
    save(fig, out, name)
    return dict(task=2, file=name, sources=source, panels=["Finite-total EPT curves for four cooperativities", "Peak ternary complex over the two-affinity grid at alpha=10", "Two 100-conformer distance histograms"],
                caveats=["Hook effect is conditional, not a universal clinical response", "Numerical optimum differs from the supplied approximate formula", "ETKDG/MMFF sample counts are not Boltzmann weights"])


def task3(repo: Path, out: Path) -> dict:
    base = "projects/task03_covalent_kinetics/data_task3/"
    source = [base+n for n in ["kobs_concentration.csv", "washout_trajectories.csv", "warhead_results.csv"]]
    kinetics, washout, warheads = [rows(repo, path) for path in source]
    eff = vals(warheads, "kinact_over_KI_QSSA_M_inv_s")
    gsh = vals(warheads, "gsh_forward_half_life_min")
    nondominated = np.array([not np.any((eff >= e) & (gsh >= g) & ((eff > e) | (gsh > g))) for e, g in zip(eff, gsh)])
    derived = out.parent / "data_omnibus" / "task3_pareto_qssa.csv"
    derived.parent.mkdir(parents=True, exist_ok=True)
    with derived.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["id", "kinact_over_KI_QSSA_M_inv_s", "gsh_forward_half_life_min", "pareto_nondominated_qssa", "evidence"], lineterminator="\n")
        writer.writeheader()
        for row, e, g, flag in zip(warheads, eff, gsh, nondominated):
            writer.writerow(dict(id=row["id"], kinact_over_KI_QSSA_M_inv_s=e, gsh_forward_half_life_min=g,
                                 pareto_nondominated_qssa=bool(flag), evidence="Derived from uncalibrated scenario inputs; not a DILI model"))
    fig, gs = canvas(3, "Covalent kinetics and the target / GSH trade-off",
        "Evidence: eight illustrative warhead scenarios; kinetic and GSH parameters are uncalibrated. Forward-only kobs panels suppress chemical reversal.\n"
        "Recovery restores each scenario's reversibility and target turnover. The Pareto set maximizes QSSA efficiency and forward GSH half-life; it does not predict DILI.")
    ax = fig.add_subplot(gs[0])
    for row, color in zip(warheads, COLORS):
        sub = select(kinetics, id=row["id"])
        x = vals(sub, "concentration_M")*1e6
        ax.plot(x, vals(sub, "kobs_QSSA_s"), color=color, label=row["id"])
        ax.scatter(x[::6], vals(sub, "kobs_exact_s")[::6], color=color, s=10, zorder=3)
    ax.set(xscale="log", yscale="log", title="A  |  Saturating inactivation", xlabel=r"Clamped free inhibitor ($\mu$M)", ylabel=r"Observed forward rate $k_{obs}$ (s$^{-1}$)")
    ax.legend(ncol=4, loc="lower right", columnspacing=.7, handlelength=1.4)
    ax.text(.03, .96, "Lines: QSSA hyperbolae\nPoints: exact slow eigenvalue", transform=ax.transAxes, fontsize=8, va="top")
    ax = fig.add_subplot(gs[1])
    for wid, label, color, style in [("REV", "Reversible control", "#4C5861", "--"), ("W01", "W01 irreversible", COLORS[0], "-"), ("W08", "W08 reversible covalent", COLORS[7], "-")]:
        sub = [r for r in select(washout, id=wid, turnover_half_life_h=24.) if float(r["time_since_washout_h"]) >= 0]
        ax.plot(vals(sub, "time_since_washout_h"), vals(sub, "free_active_fraction")*100, color=color, ls=style, label=label)
    ax.set(title="B  |  Post-washout target recovery", xlabel="Time after washout (h)", ylabel="Free active target (%)", xlim=(0, 120), ylim=(0, 104))
    ax.legend(loc="lower right")
    ax.text(.49, .46, "2 h pulse at 1 µM\n24 h protein turnover half-life\nPerfect free-drug sink", transform=ax.transAxes, fontsize=8)
    ax = fig.add_subplot(gs[2])
    order = np.argsort(eff[nondominated])
    ax.plot(eff[nondominated][order], gsh[nondominated][order], color=INK, ls=":", lw=1.3, label="QSSA Pareto set")
    for row, e, g, color in zip(warheads, eff, gsh, COLORS):
        ax.scatter(e, g, s=65, color=color, marker="s" if float(row["krev_s"]) else "o", edgecolor="white", lw=.8, zorder=3)
        offset = (-19, -12) if row["id"] == "W04" else (5, 5)
        ax.annotate(row["id"], (e, g), xytext=offset, textcoords="offset points", fontsize=8)
    ax.set(xscale="log", yscale="log", title="C  |  Scenario Pareto frontier", xlabel=r"$k_{inact}/K_{I,QSSA}$ (M$^{-1}$ s$^{-1}$)", ylabel="Forward GSH half-life at 5 mM (min)", xlim=(200, 2e5), ylim=(1, 1600))
    ax.legend(loc="lower left")
    ax.text(.03, .22, "Larger is favored on both axes\nSquares: reversible covalent", transform=ax.transAxes, fontsize=8)
    name = "fig_task3_covalent_kinetics.png"
    save(fig, out, name)
    return dict(task=3, file=name, sources=source, derived_data=["data_omnibus/task3_pareto_qssa.csv"],
                panels=["QSSA kobs hyperbolae and exact eigenvalue points for eight scenarios", "Reversible and covalent washout at 24-hour turnover", "Recomputed QSSA-efficiency / forward GSH half-life Pareto frontier"],
                caveats=["All kinetic values are illustrative", "K_I,QSSA differs from binary K_D", "Forward GSH half-life is not a clinical safety margin or DILI prediction"])


def task4(repo: Path, out: Path) -> dict:
    base = "projects/task04_pbpk/"
    source = [base+n for n in ["single_iv.csv", "single_oral.csv", "multidose_bid_7days.csv", "steady_state_bid.csv", "results_summary.json"]]
    iv, oral, multi, steady = [rows(repo, path) for path in source[:4]]
    summary = json.loads((repo/source[4]).read_text(encoding="utf-8-sig"))
    fig, gs = canvas(4, "Whole-body exposure: route, tissue and accumulation",
        "Evidence: synthetic 7-state PBPK scenario with well-stirred IVIVE and a Poulin--Theil/logD partition approximation; this is not a Rodgers--Rowland implementation.\n"
        "Dose is 100 mg per administration. Day 7 remains below the independently solved periodic steady state; modeled exposure does not establish a clinical regimen.")
    ax = fig.add_subplot(gs[0])
    for sub, label, color in [(iv, "IV bolus", COLORS[0]), (oral, "Oral", COLORS[1])]:
        mask = vals(sub, "plasma_mg_l") > 0
        ax.plot(vals(sub, "time_h")[mask], vals(sub, "plasma_mg_l")[mask], label=label, color=color)
    ax.set(yscale="log", title="A  |  Single-dose plasma exposure", xlabel="Time (h)", ylabel="Total plasma concentration (mg/L)", xlim=(0, 48), ylim=(.01, 30))
    ax.legend(loc="upper right")
    ax.text(.03, .06, f"Modeled absolute F = {summary['oral_bioavailability_auc_ratio_pct']:.1f}%\nAUC extrapolated to infinity", transform=ax.transAxes, fontsize=8)
    ax = fig.add_subplot(gs[1])
    for tissue, color in zip(["liver", "kidney", "brain", "lung", "rest"], COLORS):
        ax.plot(vals(oral, "time_h"), vals(oral, tissue+"_mg_l"), color=color, label="Peripheral" if tissue == "rest" else tissue.title())
    ax.set(title="B  |  Distribution after oral dosing", xlabel="Time (h)", ylabel="Total tissue concentration (mg/L)", xlim=(0, 48))
    ax.legend(loc="upper right", ncol=2)
    ax = fig.add_subplot(gs[2])
    ax.plot(vals(multi, "time_h")/24, vals(multi, "plasma_mg_l"), color=COLORS[0], label="100 mg every 12 h")
    # Independent periodic solution shown on the last day, not mislabeled day-7 steady state.
    for start in [144., 156.]:
        ax.plot((vals(steady, "time_h")+start)/24, vals(steady, "plasma_mg_l"), color=COLORS[2], ls="--", label="Periodic steady state" if start == 144. else None)
    ax.set(title="C  |  Seven-day BID accumulation", xlabel="Elapsed time (days)", ylabel="Total plasma concentration (mg/L)", xlim=(0, 7), ylim=(0, 2.1))
    ax.legend(loc="upper left")
    error = summary["day7_last_interval"]["state_error_vs_ss_relative"]*100
    ax.text(.03, .72, f"Day-7 state error vs steady state\n{error:.1f}% (archived convergence metric)", transform=ax.transAxes, fontsize=8)
    name = "fig_task4_pbpk_pharmacokinetics.png"
    save(fig, out, name)
    return dict(task=4, file=name, sources=source,
                panels=["100-mg IV and oral total plasma curves", "Five oral total tissue concentration curves", "Seven-day BID trajectory versus independently solved periodic steady state"],
                caveats=["Synthetic compound inputs", "Poulin-Theil/logD partition approximation, not Rodgers-Rowland", "Day 7 is not yet periodic steady state", "No clinical dose recommendation"])


def task5(repo: Path, out: Path) -> dict:
    base = "projects/task05_cyp_ddi/results/"
    source = [base+n for n in ["synthetic_preincubation.csv", "probe_pk_curves.csv", "conditional_risk_grid.csv"]]
    incubation, curves, risk = [rows(repo, path) for path in source]
    fig, gs = canvas(5, "Time-dependent CYP inhibition and dynamic DDI",
        "Evidence: named-drug synthetic scenarios with assumed numerical inputs; the midazolam 2 mg oral probe is modeled, not a clinical observation.\n"
        "M12-form screening cutoffs trigger further evaluation, not a safety verdict. The matrix assumes Iu = KI = 1 µM, fm = 0.9, Fg = 1, kdeg = 0.019 h⁻¹.")
    ax = fig.add_subplot(gs[0])
    for compound, color in [("Clarithromycin", COLORS[0]), ("Ritonavir", COLORS[2])]:
        for concentration, style in [(.1, ":"), (2., "--"), (10., "-")]:
            sub = select(incubation, compound=compound, concentration_uM=concentration)
            ax.plot(vals(sub, "time_h"), vals(sub, "relative_activity")*100, color=color, ls=style)
    handles = [Line2D([], [], color=COLORS[0], label="Clarithromycin"), Line2D([], [], color=COLORS[2], label="Ritonavir")]
    handles += [Line2D([], [], color=INK, ls=style, label=f"{c:g} µM") for c, style in [(.1, ":"), (2, "--"), (10, "-")]]
    ax.legend(handles=handles, ncol=2, loc="lower left", columnspacing=1.)
    ax.set(title="A  |  Constructed preincubation assay", xlabel="Preincubation time (h)", ylabel="Remaining enzyme activity (%)", xlim=(0, 4), ylim=(0, 103))
    ax = fig.add_subplot(gs[1])
    for compound, color in [("Control", "#56616A"), ("Clarithromycin", COLORS[0]), ("Ritonavir", COLORS[2]), ("Ketoconazole", COLORS[3])]:
        sub = [r for r in select(curves, compound=compound, probe="Midazolam", probe_time_h=312.) if float(r["time_after_probe_h"]) <= 48]
        ax.plot(vals(sub, "time_after_probe_h"), vals(sub, "plasma_mg_l")*1000, color=color, label=compound)
    ax.set(title="B  |  Day-14 midazolam probe", xlabel="Hours after probe administration", ylabel="Total midazolam plasma (ng/mL)", xlim=(0, 48))
    ax.legend(loc="upper right")
    ax = fig.add_subplot(gs[2])
    xx, yy, zz = grid(risk, "imax_u_over_ki", "kinact_over_KI_h_uM", "category")
    cmap = ListedColormap(["#B8DBCF", "#F0CF7E", "#C97574"])
    ax.pcolormesh(xx, yy, zz, cmap=cmap, norm=BoundaryNorm([-.5,.5,1.5,2.5], 3), shading="auto", rasterized=True)
    ax.axvline(.02, color=INK, ls="--", lw=1)
    ax.set(xscale="log", yscale="log", title="C  |  Conditional regulatory screen", xlabel=r"$I_{max,u}/K_i$ (dimensionless)", ylabel=r"$k_{inact}/K_I$ (h$^{-1}$ $\mu$M$^{-1}$)")
    handles = [Line2D([], [], marker="s", ls="", markersize=8, color=color, label=label) for color, label in zip(cmap.colors, ["Below basic cutoffs", "Further evaluation", "Model AUCR ≥ 5"])]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(.5, -.23), ncol=1, fontsize=8, labelspacing=.3)
    name = "fig_task5_cyp_ddi_mbi.png"
    save(fig, out, name)
    return dict(task=5, file=name, sources=source,
                panels=["Six noise-free preincubation curves for two MBI scenarios", "Hypothetical day-14 midazolam dynamic PK under three perpetrators and control", "Conditional M12-form screening matrix from archived assumptions"],
                caveats=["Drug names do not validate assumed numerical parameters", "Conditional basic screen is not an FDA clinical risk classification", "Screening ratios are not themselves AUCR", "No inference of clinical safety below cutoffs"])


def generate(repo: Path, out: Path) -> list[dict]:
    """Write five 300-DPI panels and return auditable source/caveat metadata."""
    repo, out = Path(repo).resolve(), Path(out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    with plt.rc_context(STYLE):
        return [fn(repo, out) for fn in (task1, task2, task3, task4, task5)]


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = generate(args.repo, args.out)
    (args.out.parent / "figures_1_5_metadata.json").write_text(json.dumps(result, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({"figures": len(result), "out": str(args.out)}, ensure_ascii=False))
