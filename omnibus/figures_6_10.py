"""Source-backed omnibus figures for Tasks 6--10.

The project artifacts are read-only inputs. Small plotting-derived arrays are
written next to the figure directory under ``data_omnibus`` for inspection.
No missing measurements are imputed and no clinical calibration is inferred.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm, LogNorm
import numpy as np
import pandas as pd


COLORS = ["#087F8C", "#D88532", "#6750A4", "#3267A3", "#BE4963"]
INK = "#17324D"
STYLE = {
    "font.family": "DejaVu Sans", "font.size": 10,
    "axes.titlesize": 11.5, "axes.titleweight": "bold",
    "axes.labelsize": 10, "axes.edgecolor": "#B8C4CE",
    "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": "#41566A", "ytick.color": "#41566A",
    "axes.spines.top": False, "axes.spines.right": False,
    "grid.color": "#DDE5EB", "grid.alpha": .7,
    "legend.frameon": False, "legend.fontsize": 8,
    "lines.linewidth": 2.1, "savefig.facecolor": "white",
}


def _json(path: Path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def _csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path)


def _title(ax, letter: str, title: str):
    ax.set_title(f"{letter}  {title}", loc="left", pad=11)
    ax.set_axisbelow(True)


def _finish(fig, out: Path, name: str, title: str, subtitle: str, note: str):
    fig.suptitle(title, x=.065, y=.985, ha="left", fontsize=19,
                 fontweight="bold", color=INK)
    fig.text(.065, .985-.44/fig.get_figheight(), subtitle,
             fontsize=10.2, color="#526B80", ha="left")
    fig.text(.065, .021, note, fontsize=8.5, color="#526B80", va="bottom")
    fig.savefig(out / name, dpi=300, facecolor="white")
    plt.close(fig)


def _write_derived(out: Path, name: str, rows: pd.DataFrame) -> str:
    folder = out.parent / "data_omnibus"
    folder.mkdir(parents=True, exist_ok=True)
    path = folder / name
    rows.to_csv(path, index=False, lineterminator="\n", float_format="%.12g")
    return "data_omnibus/" + name


def _meta(task, name, paths, repo, panels, caveats, derived=()):
    return {"task": task, "file": name,
            "sources": [p.relative_to(repo).as_posix() for p in paths],
            "panels": panels, "caveats": caveats, "derived_files": list(derived)}


def _task6(repo: Path, out: Path) -> dict:
    p = repo / "projects/task06_asd"
    sources = [p / "inputs/parameters.json", p / "results/phase_boundaries.csv",
               p / "results/concentration_mass_timeseries.csv", p / "results/summary.json"]
    cfg, phase, curves, summary = _json(sources[0]), _csv(sources[1]), _csv(sources[2]), _json(sources[3])
    d = cfg["drug"]
    w = np.linspace(.001, .999, 301)
    derived = []
    fig, axs = plt.subplots(2, 2, figsize=(13.5, 10.1))
    fig.subplots_adjust(left=.075, right=.965, bottom=.12, top=.865,
                        wspace=.28, hspace=.49)
    for color, polymer in zip(COLORS, cfg["polymers"]):
        phi = (w / d["density_g_cm3"]) / (w / d["density_g_cm3"] + (1-w) / polymer["density_g_cm3"])
        delta = np.asarray(d["hansen_MPa_half"]) - polymer["hansen_MPa_half"]
        chi = d["molar_volume_cm3_mol"] * np.dot(delta**2, [1, .25, .25]) / (8.31446261815324 * 298.15)
        g = phi * np.log(phi) / d["N_segments"] + (1-phi) * np.log1p(-phi) / polymer["N_segments"] + chi * phi * (1-phi)
        k = d["density_g_cm3"] * d["Tg_K"] / (polymer["density_g_cm3"] * polymer["Tg_K"])
        tg = (w*d["Tg_K"] + k*(1-w)*polymer["Tg_K"]) / (w+k*(1-w))
        axs[0, 0].plot(w, g, color=color, label=polymer["name"])
        axs[1, 0].plot(w, tg-273.15, color=color, label=polymer["name"])
        rows = phase[(phase.polymer == polymer["name"]) & phase.two_phase]
        if len(rows):
            for side in ("low", "high"):
                axs[0, 1].plot(rows[f"binodal_{side}_w"], rows.temperature_K-273.15,
                               color=color, label=polymer["name"] if side == "low" else None)
                axs[0, 1].plot(rows[f"spinodal_{side}_w"], rows.temperature_K-273.15,
                               color=color, ls="--", lw=1.35)
        derived.extend({"polymer": polymer["name"], "drug_mass_fraction": wi,
                        "drug_volume_fraction": pi, "chi_25C": chi,
                        "mixing_free_energy_over_RT_per_site": gi,
                        "dry_Tg_C": ti-273.15}
                       for wi, pi, gi, ti in zip(w, phi, g, tg))
    _title(axs[0, 0], "A", "Mixing free energy at 25 °C")
    axs[0, 0].set(xlabel="Drug mass fraction", ylabel=r"$\Delta g_{mix}/RT$ per lattice site")
    axs[0, 0].axhline(0, color="#8996A1", lw=.8)
    axs[0, 0].set_ylim(-.32, .025)
    axs[0, 0].legend(ncol=2, loc="lower left")
    _title(axs[0, 1], "B", "Amorphous phase separation envelopes")
    axs[0, 1].set(xlabel="Drug mass fraction", ylabel="Temperature (°C)", xlim=(0, 1), ylim=(0, 150))
    axs[0, 1].legend(loc="upper right")
    axs[0, 1].text(.03, .96, "Solid: binodal | Dashed: spinodal\nHPMC-AS: one phase over 0–150 °C",
                   transform=axs[0, 1].transAxes, fontsize=8, va="top",
                   bbox={"facecolor": "white", "edgecolor": "none", "alpha": .9})
    _title(axs[1, 0], "C", "Gordon–Taylor dry glass transition")
    axs[1, 0].set(xlabel="Drug mass fraction", ylabel=r"$T_{g,mix}$ (°C)")
    axs[1, 0].legend(ncol=2)
    for color, formulation, label in zip(["#8594A2", COLORS[1], COLORS[0]],
                                        ["crystalline", "amorphous", "ASD"],
                                        ["Crystalline", "Amorphous", "20% HPMC-AS ASD"]):
        rows = curves[(curves.formulation == formulation) & ~curves.absorption_enabled]
        axs[1, 1].plot(rows.time_h, rows.free_concentration_mg_ml / d["crystalline_solubility_mg_ml"], color=color, label=label)
    axs[1, 1].axhline(1, color="#8996A1", ls=":", lw=1)
    _title(axs[1, 1], "D", "Closed-vessel supersaturation")
    axs[1, 1].set(xlabel="Time (h)", ylabel=r"$C/C_{crystalline}$", xlim=(0, 6), ylim=(0, None))
    axs[1, 1].legend(loc="center right", bbox_to_anchor=(.99, .40))
    axs[1, 1].text(.45, .75, f"ASD concentration AUC / crystal = {summary['ASD20_AUC_ratio']:.2f}\n"
                   f"Time above 1.05× = {summary['ASD20_closed_vessel']['supersaturation_over_1p05_duration_h']:.2f} h",
                   transform=axs[1, 1].transAxes, va="top", fontsize=8.5)
    for ax in axs.flat:
        ax.grid(alpha=.45)
    name = "fig_task6_asd_supersaturation.png"
    file = _write_derived(out, "task6_thermodynamic_plot_grid.csv", pd.DataFrame(derived))
    _finish(fig, out, name, "06  |  Formulation thermodynamics",
            "From miscibility to the supersaturation trajectory • 4 carrier scenarios • archived mass-balanced ODE results",
            "Assumed API / Hansen / kinetic parameters; named carriers are not calibrated products.\n"
            "The medium is an intestinal-style scenario, not a validated FaSSIF experiment; concentration AUC is not systemic exposure.")
    return _meta(6, name, sources, repo,
                 ["FH mixing free energy", "archived binodal/spinodal curves", "dry Gordon–Taylor Tg", "closed-vessel dissolution"],
                 ["Hypothetical API and carrier parameterization; no DSC or dissolution calibration.",
                  "HPMC-AS has no phase separation envelope in the computed range; absent boundaries are not zeros.",
                  "20% ASD is a prespecified candidate, not a demonstrated optimum."], [file])


def _task7(repo: Path, out: Path) -> dict:
    p = repo / "projects/task07_qsp"
    sources = [p / "trajectories.csv", p / "dose_synergy_matrix.csv", p / "kaplan_meier.csv",
               p / "progression_endpoints.csv", p / "results_summary.json"]
    trajectories, doses, km, endpoints = [_csv(s) for s in sources[:4]]
    summary = _json(sources[4])
    n = summary["n_unique_virtual_people"]
    grouped = trajectories.groupby(["arm", "time_day"])["total_mass_mg"].agg(
        median="median", q10=lambda x: x.quantile(.1), q90=lambda x: x.quantile(.9)).reset_index()
    if trajectories.groupby("arm").patient.nunique().ne(n).any():
        raise ValueError("Task7 unique-person count is inconsistent")
    fig, axs = plt.subplots(1, 3, figsize=(16.2, 6.0))
    fig.subplots_adjust(left=.06, right=.975, top=.80, bottom=.24, wspace=.36)
    arms = ["Vehicle", "Targeted", "Anti-PD-1", "Combination"]
    colors = ["#8594A2", COLORS[1], COLORS[2], COLORS[0]]
    for arm, color in zip(arms, colors):
        q = grouped[grouped.arm == arm]
        axs[0].plot(q.time_day, q["median"], color=color, label=arm)
        axs[0].fill_between(q.time_day, q.q10, q.q90, color=color, alpha=.12)
        rows = km[km.arm == arm].sort_values("time_day")
        time = rows.time_day.to_numpy()
        survival = rows.survival.to_numpy()
        if time[-1] < summary["horizon_day"]:
            time = np.r_[time, summary["horizon_day"]]
            survival = np.r_[survival, survival[-1]]
        axs[2].step(time, survival, where="post", color=color, label=arm)
        censored = endpoints[(endpoints.arm == arm) & (endpoints.event == 0)]
        if len(censored):
            ct = censored.duration_day.to_numpy()
            cs = survival[np.maximum(0, np.searchsorted(time, ct, side="right")-1)]
            axs[2].plot(ct, cs, "+", color=color, ms=7)
    _title(axs[0], "A", f"Four arms · {n} matched virtual people")
    axs[0].set(xlabel="Time (day)", ylabel="Tumor mass (mg)", xlim=(0, 60), ylim=(0, None))
    axs[0].text(.025, .97, "Median and 10–90% range", transform=axs[0].transAxes, va="top", fontsize=8)
    axs[0].legend(loc="upper left", bbox_to_anchor=(0, -.19), ncol=2)
    matrix = doses.pivot(index="mab_dose_mg_kg", columns="small_dose_mg", values="excess_over_bliss")
    if matrix.isna().any().any():
        raise ValueError("Task7 Bliss matrix has missing grid cells")
    limit = max(float(np.abs(matrix.to_numpy()).max()), 1e-6)
    mesh = axs[1].imshow(matrix, origin="lower", aspect="auto", cmap="RdBu_r",
                         norm=TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit))
    axs[1].set_xticks(range(len(matrix.columns)), [f"{x:g}" for x in matrix.columns])
    axs[1].set_yticks(range(len(matrix.index)), [f"{x:g}" for x in matrix.index])
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            v = float(matrix.iloc[i, j])
            axs[1].text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=8,
                         color="white" if abs(v) > .65*limit else INK)
    _title(axs[1], "B", "Day-60 nominal Bliss interaction")
    axs[1].set(xlabel="Small-molecule dose (mg)", ylabel="Antibody dose (mg/kg)")
    cb = fig.colorbar(mesh, ax=axs[1], fraction=.044, pad=.025)
    cb.set_label("Excess over Bliss", fontsize=9)
    _title(axs[2], "C", "Time to model mass threshold")
    axs[2].set(xlabel="Time (day)", ylabel="Probability without threshold crossing",
                xlim=(0, 60), ylim=(-.025, 1.025))
    axs[2].text(.97, .94, "Event: mass > 1.2 × baseline\n+ : administrative censoring",
                 transform=axs[2].transAxes, va="top", ha="right", fontsize=8,
                 bbox={"facecolor": "white", "edgecolor": "none", "alpha": .85})
    axs[0].grid(alpha=.4)
    axs[2].grid(alpha=.4)
    name = "fig_task7_qsp_immuno_oncology.png"
    file = _write_derived(out, "task7_tumor_quantiles.csv", grouped)
    _finish(fig, out, name, "07  |  Tumor–immune systems pharmacology",
            "Simeoni-style tumor kinetics + immune recruitment • counterfactual simulation with the same 50 people in every arm",
            "Synthetic, uncalibrated population. Bliss is conditional on these response models.\n"
            "Kaplan–Meier describes a tumor-mass threshold proxy; it is not RECIST, clinical PFS, or evidence of patient benefit.")
    return _meta(7, name, sources, repo,
                 ["matched 50-person tumor trajectories", "nominal 6×6 Bliss surface", "mass-threshold Kaplan–Meier with censor marks"],
                 ["Shading is the population 10–90% interval, not a confidence interval.",
                  "The same synthetic people are reused in four arms; 200 trajectories do not mean 200 unique people.",
                  "The mass endpoint is not clinical PFS."], [file])


def _task8(repo: Path, out: Path) -> dict:
    p = repo / "projects/task08_adc"
    sources = [p / "config.json", p / "data/dar_timecourses.csv", p / "data/hic_chromatogram.csv",
               p / "data/cell_trafficking.csv", p / "data/radial_timecourses.csv"]
    cfg = _json(sources[0])
    dar, hic, cell, radial = [_csv(s) for s in sources[1:]]
    endpoint = dar.loc[dar.groupby("payload_equivalents").time_h.idxmax()]
    base = endpoint[endpoint.payload_equivalents == cfg["downstream_payload_equivalents"]]
    if len(base) != 1:
        raise ValueError("Task8 downstream DAR scenario is not uniquely identified")
    mean = float(base.mean_dar.iloc[0])
    lysosome = cell.lysosome_adc_molecules_per_cell
    flux = cfg["cleavage_vmax_adc_per_cell_h"] * lysosome / (cfg["cleavage_km_adc_per_cell"]+lysosome) * mean
    derived = cell[["permeability_scenario", "time_h"]].copy()
    derived["cleaved_payload_molecules_per_cell_h"] = flux
    fig = plt.figure(figsize=(15.7, 10.4))
    gs = fig.add_gridspec(2, 6, left=.065, right=.955, bottom=.115, top=.865,
                         wspace=.85, hspace=.50)
    axes = [fig.add_subplot(gs[0, i:i+2]) for i in (0, 2, 4)]
    axes += [fig.add_subplot(gs[1, 0:3]), fig.add_subplot(gs[1, 3:6])]
    width = .19
    x = np.arange(9)
    for j, (_, row) in enumerate(endpoint.iterrows()):
        values = np.array([row[f"dar_{i}_fraction"] for i in x])
        if not np.isclose(values.sum(), 1, atol=1e-7):
            raise ValueError("Task8 endpoint DAR fractions are not normalized")
        axes[0].bar(x+(j-1.5)*width, values, width, color=COLORS[j],
                    label=f"{row.payload_equivalents:g} eq · mean {row.mean_dar:.2f}")
    _title(axes[0], "A", "Finite-linker conjugation at 8 h")
    axes[0].set(xlabel="Drug-to-antibody ratio", ylabel="Species fraction", xticks=x)
    axes[0].legend(fontsize=7.5)
    dar_colors = plt.get_cmap("viridis")(np.linspace(.08, .9, 9))
    for j in x:
        axes[1].fill_between(hic.time_min, 0, hic[f"dar_{j}_signal"],
                             color=dar_colors[j], alpha=.28, linewidth=0)
        axes[1].plot(hic.time_min, hic[f"dar_{j}_signal"], color=dar_colors[j], lw=1,
                      label=f"DAR {j}")
    axes[1].plot(hic.time_min, hic.total_signal, color=INK, lw=2, label="Mixture")
    _title(axes[1], "B", "Synthetic HIC analytical profile")
    axes[1].set(xlabel="Retention time (min)", ylabel="Signal (normalized area/min)", xlim=(0, 28))
    axes[1].legend(ncol=2, fontsize=7, loc="upper right")
    for scenario, color, style in [("High permeability", COLORS[0], "-"), ("Low permeability", COLORS[1], "--")]:
        rows = derived[derived.permeability_scenario == scenario]
        axes[2].plot(rows.time_h, rows.cleaved_payload_molecules_per_cell_h/1e3,
                      color=color, ls=style, label=scenario)
    _title(axes[2], "C", "Lysosomal payload generation")
    axes[2].set(xlabel="Time (h)", ylabel=r"Cleavage flux ($10^3$ molecules cell$^{-1}$ h$^{-1}$)", xlim=(0, 72))
    axes[2].legend(loc="lower right", fontsize=8)
    axes[2].text(.20, .70, "Curves overlap: cleavage is upstream\nof the permeability-dependent export step.",
                 transform=axes[2].transAxes, va="top", fontsize=8)
    maximum = float(radial.extracellular_nm.max())
    for ax, letter, scenario in zip(axes[3:], ["D", "E"], ["High permeability", "Low permeability"]):
        rows = radial[radial.scenario == scenario]
        z = rows.pivot(index="time_h", columns="radius_um", values="extracellular_nm")
        if z.isna().any().any():
            raise ValueError("Task8 radial mesh has missing cells")
        mesh = ax.pcolormesh(z.columns.to_numpy(), z.index.to_numpy(), z.to_numpy(),
                             shading="auto", cmap="magma", vmin=0, vmax=maximum)
        ax.axvline(cfg["core_radius_um"], color="white", lw=1, ls="--")
        ax.set(xlabel="Radius (µm)", ylabel="Time (h)", xlim=(0, 200))
        _title(ax, letter, f"Bystander transport · {scenario.lower()}")
        cb = fig.colorbar(mesh, ax=ax, fraction=.036, pad=.022)
        cb.set_label("Extracellular payload (nM)", fontsize=9)
    name = "fig_task8_adc_multiscale.png"
    file = _write_derived(out, "task8_lysosomal_payload_flux.csv", derived)
    _finish(fig, out, name, "08  |  ADC engineering across scales",
            "Conjugation → HIC → trafficking → conservative spherical diffusion • shared concentration scale across transport panels",
            "Mechanistic scenarios with assumed parameters. The lysosomal rate is effective cleavage, not proof of Cathepsin-B exclusivity.\n"
            "Dashed line: 20 µm source core. Concentration contours are not validated bystander killing radii; tissue geometry is idealized.")
    return _meta(8, name, sources, repo,
                 ["four endpoint DAR distributions", "DAR-resolved HIC", "model-derived lysosomal flux", "high/low permeability radial contours"],
                 ["Finite-linker mass conservation is retained; downstream mixture is not a purified DAR4 fraction.",
                  "HIC peaks and effective lysosomal rates are hypothetical.",
                  "Radial extracellular concentration is not a cytotoxicity contour."], [file])


def _task9(repo: Path, out: Path) -> dict:
    p = repo / "projects/task09_cryoem_allostery/results"
    sources = [p / "conformer_latent_trajectory.csv", p / "msm.json", p / "anm_prs.npz", p / "summary.json"]
    latent, msm = _csv(sources[0]), _json(sources[1])
    summary = _json(sources[3])
    with np.load(sources[2], allow_pickle=False) as npz:
        prs = npz["prs_normalized"].copy()
        residue = npz["residue_ids"].copy()
    fig, axs = plt.subplots(1, 3, figsize=(16.2, 6.0))
    fig.subplots_adjust(left=.06, right=.952, bottom=.245, top=.80, wspace=.37)
    states = np.sort(latent.state.unique())
    for state, color in zip(states, plt.get_cmap("viridis")(np.linspace(.05, .9, len(states)))):
        q = latent[latent.state == state]
        axs[0].scatter(q.z1, q.z2, s=15, alpha=.62, color=color,
                        linewidths=0, label=f"State {state} (n={len(q)})")
    _title(axs[0], "A", "Synthetic conformation manifold")
    axs[0].set(xlabel=r"Kernel-PCA $z_1$ (arbitrary units)", ylabel=r"Kernel-PCA $z_2$ (arbitrary units)")
    axs[0].legend(loc="upper left", fontsize=7.6)
    energy = np.asarray(msm["estimated_free_energies_kcal_mol"])
    ci = np.asarray(msm["bootstrap_free_energy_95_interval"])
    pop = np.asarray(msm["estimated_stationary"])
    x = np.arange(len(energy))
    axs[1].errorbar(x, energy, yerr=np.vstack((energy-ci[:, 0], ci[:, 1]-energy)),
                    fmt="o", color=COLORS[2], capsize=5, ms=7, elinewidth=1.8)
    for i, value in enumerate(energy):
        axs[1].text(i, ci[i, 1]+.12, f"{value:.2f}", ha="center", fontsize=8.5)
    axs[1].set_xticks(x, [f"S{i}\n{100*v:.1f}%" for i, v in enumerate(pop)])
    _title(axs[1], "B", "State population free energy")
    axs[1].set(xlabel="MSM state / stationary population", ylabel=r"$-RT\ln(\pi_i/\pi_0)$ (kcal/mol)",
                ylim=(-.15, float(ci.max()+.55)))
    axs[1].text(.03, .96, "95% parametric bootstrap intervals\nDiscrete states do not specify a barrier",
                 transform=axs[1].transAxes, va="top", fontsize=8)
    mesh = axs[2].imshow(prs, origin="lower", cmap="magma", aspect="equal",
                          norm=LogNorm(vmin=1e-3, vmax=1))
    ticks = np.unique(np.linspace(0, len(residue)-1, 6, dtype=int))
    axs[2].set_xticks(ticks, [str(residue[i]) for i in ticks])
    axs[2].set_yticks(ticks, [str(residue[i]) for i in ticks])
    _title(axs[2], "C", f"ANM response matrix · {len(residue)} Cα sites")
    axs[2].set(xlabel="Perturbed residue", ylabel="Responding residue")
    cb = fig.colorbar(mesh, ax=axs[2], fraction=.047, pad=.025)
    cb.set_label("Response / perturbed-site self-response", fontsize=9)
    name = "fig_task9_cryoem_allostery.png"
    _finish(fig, out, name, "09  |  Conformation, pockets and allostery",
            f"Public T4 lysozyme X-ray endpoints + {len(latent)} synthetic conformers • reversible count estimator • elastic-network response",
            "No cryo-EM particle reconstruction or molecular-dynamics clock. The MSM is synthetic; state free energy is not an activation barrier.\n"
            f"{summary['threshold_indeterminate_frames']} pocket frames remain indeterminate. "
            f"Contact-coupling path length is {summary['path_length_A']:.3f} Å; this is not evidence of directional energy transmission.")
    return _meta(9, name, sources, repo,
                 ["500 synthetic kernel-PCA coordinates", "MSM state population free energies and bootstrap intervals", "164×164 PRS matrix"],
                 ["T4 lysozyme is a structural method example, not an oncology target or experimental cryo-EM reconstruction.",
                  "Discrete population free energies do not identify transition-state barriers.",
                  "ROI-limited pocket volumes and the short inferred path do not satisfy requested large-pocket/30 Å claims."])


def _task10(repo: Path, out: Path) -> dict:
    p = repo / "projects/task10_rna_splicing/results"
    sources = [p / "bpp_WT_model_37C.csv", p / "secondary_ensemble_summary.json",
               p / "thermodynamic_cycles.csv", p / "splicing_mass_action.csv"]
    bpp = np.loadtxt(sources[0], delimiter=",")
    secondary = [v for v in _json(sources[1]) if v["variant"] == "WT_model" and v["temperature_C"] == 37][0]
    sequence = secondary["sequence"]
    if bpp.shape != (len(sequence), len(sequence)) or not np.allclose(bpp, bpp.T):
        raise ValueError("Task10 BPP shape or symmetry failed")
    if np.any(bpp.sum(axis=1) > 1+1e-7):
        raise ValueError("Task10 base pairing probability row exceeds one")
    cycles, dose = _csv(sources[2]), _csv(sources[3])
    fig, axs = plt.subplots(1, 3, figsize=(16.2, 6.2))
    fig.subplots_adjust(left=.06, right=.973, bottom=.25, top=.80, wspace=.37)
    mesh = axs[0].imshow(bpp, origin="lower", cmap="magma", vmin=0, vmax=1,
                         extent=(.5, len(sequence)+.5, .5, len(sequence)+.5))
    _title(axs[0], "A", f"ViennaRNA ensemble · {len(sequence)}-nt model")
    axs[0].set(xlabel="Engineered sequence position", ylabel="Engineered sequence position")
    axs[0].set_xticks([1, 5, 10, 15, 20, len(sequence)])
    axs[0].set_yticks([1, 5, 10, 15, 20, len(sequence)])
    cb = fig.colorbar(mesh, ax=axs[0], fraction=.046, pad=.025)
    cb.set_label("Base-pair probability", fontsize=9)
    axs[0].text(.02, -.23, f"37 °C | ΔGensemble = {secondary['ensemble_free_energy_kcal_mol']:.2f} kcal/mol\n"
                 "Two RNA strands joined by an artificial UUCG tether",
                 transform=axs[0].transAxes, fontsize=8.2)
    x = np.arange(len(cycles))
    short = ["Planar", "Polycation", "Mismatch", "Flexible"]
    for offset, field, color, label in [(-.24, "conf_kcal", COLORS[1], r"Competent-state penalty"),
                                       (0, "intrinsic_bind_kcal", COLORS[0], "Intrinsic binding"),
                                       (.24, "net_standard_bind_kcal", COLORS[2], "Full ensemble binding")]:
        axs[1].bar(x+offset, cycles[field], width=.23, color=color, label=label)
    axs[1].axhline(0, color="#7F8C98", lw=.8)
    axs[1].set_xticks(x, short, rotation=15)
    _title(axs[1], "B", "Hypothesized thermodynamic cycles")
    axs[1].set(ylabel="Standard free energy (kcal/mol)")
    axs[1].legend(loc="lower left", bbox_to_anchor=(0, -.28), fontsize=8)
    for i, ligand in enumerate(cycles["name"]):
        for site, style in [("WT scenario", "-"), ("mutant scenario", "--")]:
            q = dose[(dose.ligand == ligand) & (dose.site == site) & (dose.total_l_M > 0)]
            if len(q) == 0:
                raise ValueError(f"Missing Task10 dose scenario: {ligand}, {site}")
            axs[2].semilogx(q.total_l_M*1e6, q.inclusion_proxy_percent, color=COLORS[i],
                             ls=style, label=short[i] if site == "WT scenario" else None)
    _title(axs[2], "C", "Finite-mass U1 recruitment proxy")
    axs[2].set(xlabel="Total ligand (µM)", ylabel="100 × target U1 occupancy (%)", ylim=(0, 102))
    axs[2].legend(loc="upper left", fontsize=8)
    axs[2].text(.97, .06, "Solid: WT scenario\nDashed: mutant scenario",
                 transform=axs[2].transAxes, ha="right", fontsize=8)
    axs[2].grid(alpha=.4)
    name = "fig_task10_rna_targeted_cadd.png"
    _finish(fig, out, name, "10  |  RNA recognition and U1 recruitment",
            "Executed partition functions + explicit eight-state mass action • structural inspiration: public SMN-C5 RNA complexes",
            "BPP describes an engineered secondary-structure model, not a native spliceosome. Energy inputs are hypotheses, not measured ligand affinities.\n"
            "U1 occupancy is an inclusion proxy only; these curves cannot predict clinical exon inclusion, dose, or efficacy.")
    return _meta(10, name, sources, repo,
                 ["ViennaRNA WT-model BPP at 37 °C", "four hypothesized free-energy cycles", "WT/mutant eight-state U1 occupancy curves"],
                 ["The engineered hairpin replaces pseudouridine with uridine and tethers the experimental two-strand motif.",
                  "Ligand-family labels are hypothetical scenarios, not actual risdiplam/aminoglycoside measurements.",
                  "No mapping from equilibrium U1 occupancy to clinical exon inclusion is established."])


def generate(repo: Path, out: Path) -> list[dict]:
    """Render five figures from curated artifacts without modifying projects."""
    repo, out = Path(repo).resolve(), Path(out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    with plt.rc_context(STYLE):
        return [fn(repo, out) for fn in (_task6, _task7, _task8, _task9, _task10)]


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    metadata = generate(args.repo, args.out)
    print(json.dumps(metadata, ensure_ascii=False, indent=2))
