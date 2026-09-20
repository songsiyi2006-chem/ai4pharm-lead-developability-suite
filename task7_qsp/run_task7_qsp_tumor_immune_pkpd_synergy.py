"""Offline, reproducible synthetic QSP experiment; no clinical efficacy calibration.

Python >=3.10; NumPy, SciPy, Matplotlib. Times are days, masses/doses mg,
concentrations mg/L, clearance L/day, CTL an arbitrary normalized density.
The same virtual people receive every arm: paired inference is mandatory.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import platform
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.integrate import solve_ivp
from scipy.optimize import brentq

DEFAULT = {
    "seed": 20260920, "n_patients": 50, "horizon_day": 60.0,
    "output_step_day": 0.25, "body_weight_kg": 70.0,
    "initial_tumor_mg": 100.0, "doubling_time_day": 12.0,
    "doubling_time_cv": 0.25, "linear_growth_mg_day": 8.0,
    "transition_psi": 20.0, "transit_tau_day": 1.5,
    "small_dose_mg": 50.0, "small_interval_day": 1.0,
    "oral_bioavailability": 0.7, "absorption_day_inv": 12.0,
    "small_clearance_l_day": 15.0, "clearance_cv": 0.30,
    "small_volume_l": 40.0, "drug_kill_l_mg_day": 0.009,
    "mab_dose_mg_kg": 10.0, "mab_interval_day": 14.0,
    "infusion_duration_day": 1.0 / 48.0,
    "mab_clearance_l_day": 0.2, "mab_volume_l": 3.0,
    "mab_kd_mg_l": 5.0, "baseline_ctl": 1.0, "baseline_ctl_cv": 0.40,
    "antigen_stim_ctl_day": 0.9, "antigen_halfmax_mg_day": 4.0,
    "ctl_loss_day_inv": 0.08, "exhaustion_day_inv": 0.35,
    "pdl1_relative": 1.0, "immune_kill_ctl_day_inv": 0.015,
    "progression_mass_ratio": 1.20, "rtol": 2e-7, "atol": 1e-9,
    "max_step_day": 0.20, "permutation_replicates": 4095,
    "bootstrap_replicates": 4000,
    "matrix_small_doses_mg": [0, 10, 25, 50, 100, 200],
    "matrix_mab_doses_mg_kg": [0, 0.1, 0.3, 1, 3, 10],
}
ARMS = ["Vehicle", "Targeted", "Anti-PD-1", "Combination"]
COLORS = ["#697583", "#3174A6", "#D48A29", "#AA3865"]
SOURCES = [
    {"id": "simeoni2004", "doi": "10.1158/0008-5472.CAN-03-2524",
     "url": "https://pubmed.ncbi.nlm.nih.gov/14871843/",
     "role": "Preclinical growth/transit architecture; no parameter values imported."},
    {"id": "recist2009", "doi": "10.1016/j.ejca.2008.10.026",
     "url": "https://recist.eortc.org/recist-1-1/",
     "role": "Endpoint distinction: diameter sum, nadir, absolute 5 mm and new lesions."},
    {"id": "loewe2018", "doi": "10.3389/fphar.2018.00031",
     "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC5808155/",
     "role": "Dose equivalence and Loewe consistency constraints."},
    {"id": "bliss2019", "doi": "10.1371/journal.pone.0224137",
     "url": "https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0224137",
     "role": "Bliss reference and uncertainty; positive score alone is not clinical proof."},
]


def validate_config(cfg):
    if set(cfg) != set(DEFAULT):
        raise ValueError(f"Configuration keys differ: {set(cfg) ^ set(DEFAULT)}")
    for key, value in cfg.items():
        values = value if isinstance(value, list) else [value]
        if not values or any(not isinstance(v, (int, float)) or not np.isfinite(v)
                             or v < 0 for v in values):
            raise ValueError(f"Nonfinite or negative {key}")
    positive = ["n_patients", "horizon_day", "output_step_day", "body_weight_kg",
                "initial_tumor_mg", "doubling_time_day", "linear_growth_mg_day",
                "transition_psi", "transit_tau_day", "small_interval_day",
                "absorption_day_inv", "small_clearance_l_day", "small_volume_l",
                "mab_interval_day", "infusion_duration_day", "mab_clearance_l_day",
                "mab_volume_l", "mab_kd_mg_l", "antigen_halfmax_mg_day", "rtol",
                "atol", "max_step_day", "permutation_replicates", "bootstrap_replicates"]
    if any(cfg[k] <= 0 for k in positive):
        raise ValueError("Positive physical or numerical parameter required")
    if cfg["n_patients"] < 2 or not isinstance(cfg["n_patients"], int):
        raise ValueError("At least two integer patients required")
    if cfg["oral_bioavailability"] > 1 or cfg["infusion_duration_day"] >= cfg["mab_interval_day"]:
        raise ValueError("Invalid bioavailability or overlapping infusion")
    if cfg["progression_mass_ratio"] <= 1:
        raise ValueError("Progression ratio must exceed baseline")
    for key in ["seed", "permutation_replicates", "bootstrap_replicates"]:
        if not isinstance(cfg[key], int):
            raise ValueError(f"{key} must be an integer")
    for k in ["matrix_small_doses_mg", "matrix_mab_doses_mg_kg"]:
        if cfg[k][0] != 0 or np.any(np.diff(cfg[k]) <= 0):
            raise ValueError("Dose axes must increase strictly from zero")


def lognormal_mean_cv(rng, mean, cv, n):
    if mean == 0:
        return np.zeros(n)  # Degenerate zero distribution; log(0) is not evaluated.
    sigma = np.sqrt(np.log1p(cv * cv))
    return rng.lognormal(np.log(mean) - 0.5 * sigma * sigma, sigma, n)


def make_cohort(cfg, n=None, nominal=False):
    n = cfg["n_patients"] if n is None else n
    rng = np.random.default_rng(cfg["seed"])
    fields = [("clearance", "small_clearance_l_day", "clearance_cv"),
              ("doubling_time", "doubling_time_day", "doubling_time_cv"),
              ("ctl_baseline", "baseline_ctl", "baseline_ctl_cv")]
    return {name: (np.full(n, cfg[mean]) if nominal else
                   lognormal_mean_cv(rng, cfg[mean], cfg[cv], n))
            for name, mean, cv in fields}


def dose_times(interval, horizon):
    # Exclude a dose at the study end: it has no subsequent observation time.
    return np.arange(0.0, horizon - 1e-10, interval)


def oral_concentration(t, dose_mg, clearance, cfg):
    elapsed = t - dose_times(cfg["small_interval_day"], cfg["horizon_day"])
    elapsed = elapsed[elapsed >= 0]
    ke = np.atleast_1d(clearance)[:, None] / cfg["small_volume_l"]
    ka = cfg["absorption_day_inv"]
    gap = ka - ke
    same = np.abs(gap) < 1e-8
    safe_gap = np.where(same, 1.0, gap)
    kernel = ka * (np.exp(-ke * elapsed) - np.exp(-ka * elapsed)) / safe_gap
    kernel = np.where(same, ka * elapsed * np.exp(-ka * elapsed), kernel)
    return np.asarray(dose_mg) * cfg["oral_bioavailability"] / cfg["small_volume_l"] * kernel.sum(axis=1)


def antibody_concentration(t, dose_mg_kg, cfg):
    elapsed = t - dose_times(cfg["mab_interval_day"], cfg["horizon_day"])
    elapsed = elapsed[elapsed >= 0]
    dur = cfg["infusion_duration_day"]
    ke = cfg["mab_clearance_l_day"] / cfg["mab_volume_l"]
    during = np.minimum(elapsed, dur)
    post = np.maximum(0.0, elapsed - dur)
    kernel = -np.expm1(-ke * during) * np.exp(-ke * post)
    return (np.asarray(dose_mg_kg) * cfg["body_weight_kg"] /
            (dur * cfg["mab_clearance_l_day"]) * kernel.sum())


def growth(w0, lambda0, cfg):
    # Stable at huge w: lambda1 is an asymptotic mg/day rate, not a capacity.
    ratio = np.maximum(lambda0 * w0 / cfg["linear_growth_mg_day"], 1e-300)
    denominator_log = np.logaddexp(0.0, cfg["transition_psi"] * np.log(ratio)) / cfg["transition_psi"]
    return lambda0 * w0 * np.exp(-denominator_log)


def rhs(t, flat, cohort, small_dose, mab_dose, cfg):
    w0, w1, w2, w3, ctl = flat.reshape(5, -1)
    cs = oral_concentration(t, small_dose, cohort["clearance"], cfg)
    cm = antibody_concentration(t, mab_dose, cfg)
    occupancy = cm / (cfg["mab_kd_mg_l"] + cm)
    damage = cfg["drug_kill_l_mg_day"] * cs * w0
    immune = cfg["immune_kill_ctl_day_inv"] * ctl * w0
    flux = w3 / cfg["transit_tau_day"]
    exhaustion = cfg["exhaustion_day_inv"] * cfg["pdl1_relative"]
    # Basal recruitment makes pre-existing CTL and antibody monotherapy meaningful.
    basal = (cfg["ctl_loss_day_inv"] + exhaustion) * cohort["ctl_baseline"]
    stimulation = cfg["antigen_stim_ctl_day"] * flux / (cfg["antigen_halfmax_mg_day"] + flux)
    dctl = basal + stimulation - (cfg["ctl_loss_day_inv"] + exhaustion * (1 - occupancy)) * ctl
    lam0 = np.log(2) / cohort["doubling_time"]
    return np.array([growth(w0, lam0, cfg) - damage - immune,
                     damage - w1 / cfg["transit_tau_day"],
                     (w1 - w2) / cfg["transit_tau_day"],
                     (w2 - w3) / cfg["transit_tau_day"], dctl]).ravel()


def simulate(cfg, cohort, small_dose, mab_dose, *, times=None):
    n = len(cohort["clearance"])
    if times is None:
        times = np.unique(np.r_[np.arange(0, cfg["horizon_day"], cfg["output_step_day"]), cfg["horizon_day"]])
    y = np.zeros((5, n))
    y[0] = cfg["initial_tumor_mg"]
    y[4] = cohort["ctl_baseline"]
    output = np.full((len(times), 5, n), np.nan)
    output[0] = y
    oral = dose_times(cfg["small_interval_day"], cfg["horizon_day"])
    antibody = dose_times(cfg["mab_interval_day"], cfg["horizon_day"])
    cuts = np.unique(np.r_[0, oral, antibody, antibody + cfg["infusion_duration_day"], cfg["horizon_day"]])
    cuts = cuts[(cuts >= 0) & (cuts <= cfg["horizon_day"])]
    evaluations = 0
    for start, end in zip(cuts[:-1], cuts[1:]):
        sol = solve_ivp(rhs, (start, end), y.ravel(), args=(cohort, small_dose, mab_dose, cfg),
                        rtol=cfg["rtol"], atol=cfg["atol"], max_step=cfg["max_step_day"], dense_output=True)
        if not sol.success:
            raise RuntimeError(f"ODE integration failed on {start}:{end}: {sol.message}")
        indices = np.flatnonzero((times > start) & (times <= end))
        if len(indices):
            output[indices] = sol.sol(times[indices]).T.reshape(-1, 5, n)
        y = sol.y[:, -1].reshape(5, n)
        evaluations += sol.nfev
    if not np.all(np.isfinite(output)) or np.min(output) < -1e-7:
        raise RuntimeError("Invalid/nonphysical ODE output")
    return times, output, evaluations


def progression(times, mass, ratio):
    """First linearly interpolated crossing on a disclosed observation grid."""
    threshold = mass[0] * ratio
    n = mass.shape[1]
    durations = np.full(n, times[-1])
    events = np.zeros(n, dtype=bool)
    for i in range(n):
        hits = np.flatnonzero(mass[:, i] > threshold[i])
        if len(hits):
            j = hits[0]
            if j == 0:
                raise ValueError("Baseline already exceeds its own progression threshold")
            fraction = (threshold[i] - mass[j - 1, i]) / (mass[j, i] - mass[j - 1, i])
            durations[i] = times[j - 1] + fraction * (times[j] - times[j - 1])
            events[i] = True
    return durations, events


def kaplan_meier(durations, events, horizon):
    rows = [{"time_day": 0.0, "survival": 1.0, "at_risk": len(durations), "events": 0, "censored": 0}]
    survival = 1.0
    for t in np.unique(durations):
        at_risk = int(np.sum(durations >= t))
        event = int(np.sum((durations == t) & events))
        censored = int(np.sum((durations == t) & ~events))
        survival *= 1 - event / at_risk
        rows.append({"time_day": float(t), "survival": float(survival),
                     "at_risk": at_risk, "events": event, "censored": censored})
    if rows[-1]["time_day"] < horizon:
        rows.append({"time_day": horizon, "survival": survival, "at_risk": 0, "events": 0, "censored": 0})
    return rows


def paired_survival_statistics(mono, combo, cfg):
    t = np.r_[mono[0], combo[0]]
    e = np.r_[mono[1], combo[1]]
    n = len(mono[0])
    group = np.r_[np.zeros(n), np.ones(n)]
    event_times = np.unique(t[e])
    # Pooled log-rank score weights. Permuting each pair preserves its shared biology.
    weights = e.astype(float)
    for tj in event_times:
        risk = t >= tj
        weights[risk] -= np.sum(e & (t == tj)) / risk.sum()
    observed = float(weights[n:].sum())
    rng = np.random.default_rng(cfg["seed"] + 100)
    signs = rng.choice([-1, 1], (cfg["permutation_replicates"], n))
    permuted = 0.5 * (signs @ (weights[n:] - weights[:n]))
    p = (1 + np.sum(np.abs(permuted) >= abs(observed) - 1e-12)) / (len(permuted) + 1)
    delta = combo[0] - mono[0]  # Identical administrative censoring at horizon.
    bootstrap = rng.choice(delta, (cfg["bootstrap_replicates"], n), replace=True).mean(axis=1)
    result = {"design": "paired counterfactual synthetic people; not a randomized clinical trial",
              "logrank_score_combo": observed, "paired_permutation_p": float(p),
              "permutations": len(permuted), "restricted_mean_delay_day": float(delta.mean()),
              "paired_bootstrap_95pct_delay_day": np.quantile(bootstrap, [0.025, 0.975]).tolist(),
              "events_mono": int(mono[1].sum()), "events_combo": int(combo[1].sum()),
              "cox_hr_combo_vs_mono": None, "cox_cluster_robust_95pct": None,
              "cox_status": "not_estimable", "ph_assumption": "not validated; descriptive only"}
    if not np.any(mono[1]) or not np.any(combo[1]):
        result["cox_reason"] = "One or both arms have no events; finite HR not estimated."
        return result

    def score_info(beta, contributions=False):
        expbx = np.exp(beta * group)
        score, info = 0.0, 0.0
        residual = np.zeros(2 * n)
        for tj in event_times:
            risk = t >= tj
            dead = e & (t == tj)
            d = dead.sum()
            risk_sum = expbx[risk].sum()
            mean = np.dot(expbx[risk], group[risk]) / risk_sum
            score += group[dead].sum() - d * mean
            info += d * mean * (1 - mean)
            residual[dead] += group[dead] - mean
            residual[risk] -= d * expbx[risk] / risk_sum * (group[risk] - mean)
        return score, info, residual

    if score_info(-20)[0] * score_info(20)[0] >= 0:
        result["cox_reason"] = "Partial-likelihood separation or insufficient overlap."
        return result
    beta = brentq(lambda b: score_info(b)[0], -20, 20)
    _, info, residual = score_info(beta)
    cluster_scores = residual[:n] + residual[n:]
    if info <= 1e-12:
        result["cox_reason"] = "Insufficient information."
        return result
    se = math.sqrt(float(np.sum(cluster_scores ** 2) * n / (n - 1))) / info
    result.update(cox_hr_combo_vs_mono=float(np.exp(beta)),
                  cox_cluster_robust_95pct=np.exp([beta - 1.96 * se, beta + 1.96 * se]).tolist(),
                  cox_status="estimated_descriptive_cluster_robust", cox_log_hr=float(beta))
    return result


def invert_monotherapy(doses, effects, effect):
    doses, effects = np.asarray(doses), np.asarray(effects)
    if np.any(np.diff(effects) < -1e-7):
        return None, "nonmonotone_monotherapy"
    if effect <= 0:
        return None, "zero_effect_unidentifiable"
    if effect < effects[0] or effect > effects[-1]:
        return None, "outside_monotherapy_effect_support"
    # Duplicate effect plateaus have no unique inverse and are not silently extrapolated.
    exact = np.flatnonzero(np.abs(effects - effect) < 1e-12)
    if len(exact) > 1:
        return None, "nonunique_inverse"
    if len(exact) == 1:
        return float(doses[exact[0]]), "interpolated_in_support"
    j = min(len(doses) - 2, max(0, int(np.searchsorted(effects, effect)) - 1))
    if effects[j + 1] - effects[j] <= 1e-10:
        return None, "nonunique_inverse"
    return float(np.interp(effect, effects[j:j + 2], doses[j:j + 2])), "interpolated_in_support"


def loewe_index(da, db, effect, doses_a, effects_a, doses_b, effects_b):
    if da == 0 and db == 0:
        return None, "zero_dose_zero_effect"
    inverse_a, reason_a = invert_monotherapy(doses_a, effects_a, effect)
    inverse_b, reason_b = invert_monotherapy(doses_b, effects_b, effect)
    # On a monotherapy axis, only its own inverse is needed.
    if da == 0:
        return ((db / inverse_b, "axis_identity") if inverse_b else (None, reason_b))
    if db == 0:
        return ((da / inverse_a, "axis_identity") if inverse_a else (None, reason_a))
    if inverse_a is None or inverse_b is None:
        return None, f"A:{reason_a};B:{reason_b}"
    return da / inverse_a + db / inverse_b, "interpolated_in_support"


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def numeric_checks(cfg):
    checks = []
    def record(name, passed, **details):
        checks.append({"test": name, "passed": bool(passed), **details})
    huge = float(growth(np.array([1e9]), np.array([0.05]), cfg)[0])
    record("simeoni_asymptotic_linear_rate", abs(huge - cfg["linear_growth_mg_day"]) < 1e-8)
    c = make_cohort(cfg, 1, nominal=True)
    y = np.array([100, 2, 3, 4, 1.0])
    d = rhs(1, y, c, 50, 10, cfg)
    expected = growth(np.array([100]), np.log(2) / c["doubling_time"], cfg)[0] - 4 / cfg["transit_tau_day"] - cfg["immune_kill_ctl_day_inv"] * 100
    record("transit_mass_balance", abs(sum(d[:4]) - expected) < 1e-10)
    d0 = rhs(0, np.array([100, 0, 0, 0, cfg["baseline_ctl"]]), c, 0, 0, cfg)
    record("untreated_ctl_baseline_equilibrium", abs(d0[4]) < 1e-12)
    at_start = antibody_concentration(0, 10, cfg)
    at_end = antibody_concentration(cfg["infusion_duration_day"], 10, cfg)
    record("finite_infusion_zero_initial_positive_end", at_start == 0 and at_end > 0)
    equal_ke = cfg["absorption_day_inv"] * cfg["small_volume_l"]
    record("equal_absorption_elimination_limit_finite", np.isfinite(oral_concentration(0.1, 50, [equal_ke], cfg)).all())
    ci, status = loewe_index(0.25, 0.25, 0.5, [0, 1], [0, 1], [0, 1], [0, 1])
    record("loewe_dose_equivalence_identity", ci == 1.0)
    missing, _ = loewe_index(1, 1, 0.9, [0, 1], [0, 0.8], [0, 1], [0, 0.8])
    record("loewe_extrapolation_prohibited", missing is None)
    t, ev = progression(np.array([0.0, 10.0]), np.array([[100.0, 100], [140.0, 110]]), 1.2)
    record("progression_crossing_and_censoring", np.allclose(t, [5, 10]) and np.array_equal(ev, [True, False]))
    km = kaplan_meier(np.array([2., 2., 4.]), np.array([True, False, True]), 4)
    record("km_tied_death_before_censoring", abs(km[1]["survival"] - 2 / 3) < 1e-12 and km[-1]["survival"] == 0)
    null = paired_survival_statistics((t, ev), (t.copy(), ev.copy()), cfg)
    record("paired_logrank_identical_arms_null", null["paired_permutation_p"] == 1.0 and null["restricted_mean_delay_day"] == 0)
    if not all(check["passed"] for check in checks):
        raise AssertionError(checks)
    return checks


def create_figures(out, cfg, times, arrays, bliss, dose_rows, endpoints, statistics):
    folder = out / "figures_task7"
    folder.mkdir()
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, ax = plt.subplots(figsize=(9.3, 5.8), layout="constrained")
    for arm, color, values in zip(ARMS, COLORS, arrays):
        mass = values[:, :4].sum(axis=1)
        mean, sem = mass.mean(axis=1), mass.std(axis=1, ddof=1) / np.sqrt(mass.shape[1])
        ax.plot(times, mean, color=color, label=arm, lw=2)
        ax.fill_between(times, mean - sem, mean + sem, color=color, alpha=0.18)
    ax.axhline(cfg["initial_tumor_mg"] * cfg["progression_mass_ratio"], color="gray", ls=":", lw=1, label="Model progression threshold")
    ax.set(xlabel="Time (day)", ylabel="Total tumor mass (mg)", title=f"Synthetic QSP cohort: mean tumor burden ± SEM ({cfg['n_patients']} paired people)")
    ax.legend(ncol=2, fontsize=9)
    fig.savefig(folder / "fig1_simeoni_tgi_monotherapy_vs_combo.png", dpi=300)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 5.1), layout="constrained")
    for arm, color, values in zip(ARMS, COLORS, arrays):
        for ax, data in zip(axes, [values[:, 4], values[:, 3] / cfg["transit_tau_day"]]):
            avg, sem = data.mean(axis=1), data.std(axis=1, ddof=1) / np.sqrt(data.shape[1])
            ax.plot(times, avg, color=color, label=arm, lw=1.9)
            ax.fill_between(times, avg - sem, avg + sem, color=color, alpha=0.15)
    axes[0].set(xlabel="Time (day)", ylabel="Active CTL (normalized density)", title="Basal + antigen-stimulated recruitment")
    axes[1].set(xlabel="Time (day)", ylabel="Delayed cell-death flux (mg/day)", title="Antigen proxy = w3 / tau; no assay calibration")
    axes[1].legend(fontsize=9)
    fig.suptitle("Synthetic immune dynamics; mean ± SEM")
    fig.savefig(folder / "fig2_ctl_immune_infiltration_dynamics.png", dpi=300)
    plt.close(fig)
    da, db = cfg["matrix_small_doses_mg"], cfg["matrix_mab_doses_mg_kg"]
    shape = (len(db), len(da))
    eob = np.array([r["excess_over_bliss"] for r in dose_rows]).reshape(shape)
    ci = np.array([np.nan if r["loewe_ci"] is None else r["loewe_ci"] for r in dose_rows]).reshape(shape)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.8), layout="constrained")
    limit = max(0.01, np.nanmax(abs(eob)))
    image = axes[0].imshow(eob, origin="lower", cmap="RdBu_r", vmin=-limit, vmax=limit, aspect="auto")
    fig.colorbar(image, ax=axes[0], shrink=0.8, label="Excess over Bliss (fraction)")
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("#D8D8D8")
    im = axes[1].imshow(ci, origin="lower", cmap=cmap, vmin=0, vmax=max(1.5, np.nanmax(ci)), aspect="auto")
    fig.colorbar(im, ax=axes[1], shrink=0.8, label="Loewe combination index")
    for ax in axes:
        ax.set_xticks(range(len(da)), da)
        ax.set_yticks(range(len(db)), db)
        ax.set(xlabel="Oral targeted dose (mg QD)", ylabel="Antibody dose (mg/kg Q2W)")
    for j in range(shape[0]):
        for i in range(shape[1]):
            axes[0].text(i, j, f"{eob[j, i]:.2f}", ha="center", va="center", fontsize=8,
                         color="white" if abs(eob[j, i]) > limit * 0.60 else "black")
            axes[1].text(i, j, "NA" if np.isnan(ci[j, i]) else f"{ci[j, i]:.2f}", ha="center", va="center", fontsize=8,
                         color="black" if np.isnan(ci[j, i]) else "white")
    axes[0].set_title(f"Day-{cfg['horizon_day']:g} nominal-person Bliss landscape")
    axes[1].set_title("Loewe inversion; gray = unsupported / undefined")
    fig.suptitle("Synthetic dose comparison; no toxicity or clinical dose optimization")
    fig.savefig(folder / "fig3_bliss_synergy_matrix_heatmap.png", dpi=300)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.7), gridspec_kw={"width_ratios": [1.45, 1]}, layout="constrained")
    for arm, color, endpoint in zip(ARMS, COLORS, endpoints):
        km = kaplan_meier(*endpoint, cfg["horizon_day"])
        median = next((row["time_day"] for row in km if row["survival"] <= 0.5), None)
        medtext = "NR" if median is None else f"{median:.1f} d"
        axes[0].step([r["time_day"] for r in km], [r["survival"] for r in km], where="post", color=color,
                     lw=2, label=f"{arm}; median {medtext}")
        if np.any(~endpoint[1]):
            axes[0].scatter([cfg["horizon_day"]], [km[-1]["survival"]], marker="+", color=color, s=60)
    axes[0].set(xlabel="Time (day)", ylabel="Probability without model progression", ylim=(-0.03, 1.04),
                title=f"{cfg['n_patients']} paired virtual people; administrative censoring at day {cfg['horizon_day']:g}")
    axes[0].legend(fontsize=8.5, loc="upper right")
    lines = ["Combination vs targeted", "", f"Paired log-rank permutation p = {statistics['Targeted']['paired_permutation_p']:.4g}"]
    for mono in ["Targeted", "Anti-PD-1"]:
        s = statistics[mono]
        hr = s["cox_hr_combo_vs_mono"]
        hrtext = "not estimable" if hr is None else f"{hr:.3f} (cluster-robust)"
        lines.extend(["", f"Vs {mono}:", f"  HR: {hrtext}", f"  Events: {s['events_mono']} vs {s['events_combo']}",
                      f"  Restricted mean delay: {s['restricted_mean_delay_day']:.2f} d"])
    lines.extend(["", f"Mass >{100 * cfg['progression_mass_ratio']:g}% baseline is a model endpoint.", "It is NOT RECIST 1.1 or clinical PFS.", "HR assumes PH; not validated here.", "All statistical uncertainty is simulation-only."])
    axes[1].axis("off")
    axes[1].text(0, 1, "\n".join(lines), ha="left", va="top", transform=axes[1].transAxes, fontsize=9.3)
    fig.savefig(folder / "fig4_virtual_cohort_kaplan_meier_pfs.png", dpi=300)
    plt.close(fig)


def run(cfg, out):
    validate_config(cfg)
    if out.exists():
        raise FileExistsError("Use a new --out directory to preserve prior evidence.")
    out.mkdir(parents=True)
    started = time.perf_counter()
    write_json(out / "parameters_used.json", {"evidence_type": "synthetic_model_assumptions", "parameters": cfg})
    write_json(out / "sources.json", {"accessed_date": "2026-09-20", "sources": SOURCES,
                                      "clinical_data": None, "experimental_data": None})
    checks = numeric_checks(cfg)
    cohort = make_cohort(cfg)
    write_csv(out / "virtual_patients.csv", [dict(patient=i + 1, evidence="synthetic_fixed_seed",
                  clearance_l_day=float(cohort["clearance"][i]), doubling_time_day=float(cohort["doubling_time"][i]),
                  ctl_baseline=float(cohort["ctl_baseline"][i])) for i in range(cfg["n_patients"])])
    arrays, endpoints, rows, km_rows, simulation_log = [], [], [], [], []
    regimes = [(0, 0), (cfg["small_dose_mg"], 0), (0, cfg["mab_dose_mg_kg"]),
               (cfg["small_dose_mg"], cfg["mab_dose_mg_kg"])]
    for arm, (da, db) in zip(ARMS, regimes):
        t0 = time.perf_counter()
        times, values, nfev = simulate(cfg, cohort, da, db)
        arrays.append(values)
        mass = values[:, :4].sum(axis=1)
        endpoint = progression(times, mass, cfg["progression_mass_ratio"])
        endpoints.append(endpoint)
        simulation_log.append(dict(arm=arm, n=cfg["n_patients"], nfev=nfev, elapsed_seconds=time.perf_counter() - t0))
        for i, day in enumerate(times):
            for p in range(cfg["n_patients"]):
                rows.append(dict(arm=arm, patient=p + 1, time_day=float(day),
                                 w0_mg=float(values[i, 0, p]), w1_mg=float(values[i, 1, p]),
                                 w2_mg=float(values[i, 2, p]), w3_mg=float(values[i, 3, p]),
                                 total_mass_mg=float(mass[i, p]), ctl=float(values[i, 4, p]),
                                 antigen_proxy_mg_day=float(values[i, 3, p] / cfg["transit_tau_day"])))
        for row in kaplan_meier(*endpoint, cfg["horizon_day"]):
            km_rows.append(dict(arm=arm, **row))
        print(f"Completed {arm}: {int(endpoint[1].sum())}/{cfg['n_patients']} model progression events", flush=True)
    write_csv(out / "trajectories.csv", rows)
    write_csv(out / "kaplan_meier.csv", km_rows)
    write_csv(out / "progression_endpoints.csv", [dict(arm=arm, patient=p + 1, duration_day=float(endpoint[0][p]),
              event=int(endpoint[1][p]), endpoint=f"mass_gt_{100 * cfg['progression_mass_ratio']:g}pct_baseline_not_RECIST")
              for arm, endpoint in zip(ARMS, endpoints) for p in range(cfg["n_patients"])])
    masses = [values[:, :4].sum(axis=1) for values in arrays]
    effects = [1 - mass / masses[0] for mass in masses]
    bliss = effects[3] - (effects[1] + effects[2] - effects[1] * effects[2])
    write_csv(out / "bliss_timecourse.csv", [dict(time_day=float(t), mean_eob=float(bliss[i].mean()),
              sem_eob=float(bliss[i].std(ddof=1) / np.sqrt(cfg["n_patients"])),
              mean_targeted_effect=float(effects[1][i].mean()), mean_mab_effect=float(effects[2][i].mean()),
              mean_combo_effect=float(effects[3][i].mean())) for i, t in enumerate(times)])
    statistics = {ARMS[k]: paired_survival_statistics(endpoints[k], endpoints[3], cfg) for k in [1, 2]}
    write_json(out / "survival_statistics.json", statistics)

    # A dose landscape is defined for one nominal parameter set, separately from cohort uncertainty.
    da = np.array(cfg["matrix_small_doses_mg"], dtype=float)
    db = np.array(cfg["matrix_mab_doses_mg_kg"], dtype=float)
    aa, bb = np.meshgrid(da, db)
    inv_a = np.unique(np.r_[da, np.geomspace(0.01, 800, 33)])
    inv_b = np.unique(np.r_[db, np.geomspace(0.001, 80, 33)])
    doses_a = np.r_[aa.ravel(), inv_a, np.zeros(len(inv_b))]
    doses_b = np.r_[bb.ravel(), np.zeros(len(inv_a)), inv_b]
    nominal = make_cohort(cfg, len(doses_a), nominal=True)
    _, grid_values, grid_nfev = simulate(cfg, nominal, doses_a, doses_b)
    grid_mass = grid_values[-1, :4].sum(axis=0)
    nominal_effect = 1 - grid_mass / grid_mass[0]
    number_matrix = aa.size
    mono_a = nominal_effect[number_matrix:number_matrix + len(inv_a)]
    mono_b = nominal_effect[number_matrix + len(inv_a):]
    dose_rows = []
    for k, (a, b, combined) in enumerate(zip(aa.ravel(), bb.ravel(), nominal_effect[:number_matrix])):
        ea, eb = np.interp(a, inv_a, mono_a), np.interp(b, inv_b, mono_b)
        ci, status = loewe_index(a, b, combined, inv_a, mono_a, inv_b, mono_b)
        dose_rows.append(dict(small_dose_mg=float(a), mab_dose_mg_kg=float(b), day=cfg["horizon_day"],
                             targeted_effect=float(ea), mab_effect=float(eb), combo_effect=float(combined),
                             excess_over_bliss=float(combined - ea - eb + ea * eb), loewe_ci=ci, loewe_status=status))
    write_csv(out / "dose_synergy_matrix.csv", dose_rows)
    write_csv(out / "loewe_monotherapy_support.csv", [dict(agent=agent, dose=float(dose), endpoint_day=cfg["horizon_day"], endpoint_effect=float(effect))
              for agent, ds, es in [("targeted_mg_QD", inv_a, mono_a), ("mab_mg_kg_Q2W", inv_b, mono_b)]
              for dose, effect in zip(ds, es)])
    # PK is retained explicitly, including the initial half-hour infusion profile.
    pk_times = np.unique(np.r_[times, np.linspace(0, cfg["infusion_duration_day"], 9)])
    pk_rows = []
    for t in pk_times:
        cs = oral_concentration(t, cfg["small_dose_mg"], cohort["clearance"], cfg)
        cm = float(antibody_concentration(t, cfg["mab_dose_mg_kg"], cfg))
        pk_rows.append(dict(time_day=float(t), small_mean_mg_l=float(cs.mean()), small_min_mg_l=float(cs.min()),
                            small_max_mg_l=float(cs.max()), mab_mg_l=cm, receptor_occupancy=cm / (cfg["mab_kd_mg_l"] + cm)))
    write_csv(out / "pk_exposure.csv", pk_rows)

    # Verify solver and event-grid sensitivity on a single paired nominal person.
    nominal_one = make_cohort(cfg, 1, nominal=True)
    tighter = dict(cfg, rtol=cfg["rtol"] / 10, atol=cfg["atol"] / 10, max_step_day=cfg["max_step_day"] / 2,
                   output_step_day=cfg["output_step_day"] / 2)
    sensitivity_rows = []
    for arm, (a, b) in zip(ARMS, regimes):
        t, coarse, _ = simulate(cfg, nominal_one, a, b)
        tf, fine, _ = simulate(tighter, nominal_one, a, b)
        mc, mf = coarse[:, :4].sum(axis=1), fine[:, :4].sum(axis=1)
        relative_error = abs(mc[-1, 0] - mf[-1, 0]) / max(mf[-1, 0], 1e-10)
        ec, ef = progression(t, mc, cfg["progression_mass_ratio"]), progression(tf, mf, cfg["progression_mass_ratio"])
        sensitivity_rows.append(dict(arm=arm, endpoint_day=cfg["horizon_day"], relative_endpoint_mass_error=float(relative_error),
                                     coarse_progression_day=float(ec[0][0]), fine_progression_day=float(ef[0][0]),
                                     event_time_difference_day=float(abs(ec[0][0] - ef[0][0]))))
    checks.append(dict(test="nominal_solver_and_grid_refinement", passed=all(r["relative_endpoint_mass_error"] < 1e-4 and r["event_time_difference_day"] < 0.01 for r in sensitivity_rows), details=sensitivity_rows))
    if not all(c["passed"] for c in checks):
        raise AssertionError(checks)
    write_json(out / "verification.json", {"checks": checks, "scope": "numerical/software only; not biological or clinical validation"})
    create_figures(out, cfg, times, arrays, bliss, dose_rows, endpoints, statistics)
    summary = {"evidence_type": "executed_synthetic_QSP_simulation", "clinical_validation": False,
               "n_unique_virtual_people": cfg["n_patients"], "n_paired_arm_trajectories": cfg["n_patients"] * 4,
               "horizon_day": cfg["horizon_day"], "timepoints": len(times), "trajectory_rows": len(rows),
               "nominal_dose_matrix_size": [len(db), len(da)], "nominal_monotherapy_support_simulations": len(inv_a) + len(inv_b),
               "grid_nfev": grid_nfev, "endpoint_day": cfg["horizon_day"], "eob_endpoint_mean": float(bliss[-1].mean()),
               "eob_endpoint_sem": float(bliss[-1].std(ddof=1) / np.sqrt(cfg["n_patients"])),
               "eob_nominal_matrix_min": min(r["excess_over_bliss"] for r in dose_rows),
               "eob_nominal_matrix_max": max(r["excess_over_bliss"] for r in dose_rows),
               "loewe_interior_estimable": sum(r["small_dose_mg"] > 0 and r["mab_dose_mg_kg"] > 0 and r["loewe_ci"] is not None for r in dose_rows),
               "loewe_interior_total": (len(da) - 1) * (len(db) - 1),
               "arms": {arm: {"endpoint_mass_mean_mg": float(mass[-1].mean()),
                                "endpoint_mass_sem_mg": float(mass[-1].std(ddof=1) / np.sqrt(cfg["n_patients"])),
                                "progression_events": int(ep[1].sum()), "administratively_censored": int((~ep[1]).sum()),
                                "median_model_progression_day": next((r["time_day"] for r in kaplan_meier(*ep, cfg["horizon_day"]) if r["survival"] <= 0.5), None)}
                        for arm, mass, ep in zip(ARMS, masses, endpoints)},
               "survival_comparisons": statistics,
               "interpretation": "Model-dependent interaction and mass-based progression; not true clinical synergy, RECIST PFS or treatment advice."}
    write_json(out / "results_summary.json", summary)
    write_json(out / "run_log.json", {"status": "complete", "elapsed_seconds": time.perf_counter() - started,
                                      "arm_integration": simulation_log, "python": platform.python_version(),
                                      "numpy": np.__version__, "scipy": scipy.__version__, "matplotlib": matplotlib.__version__})
    script_hash = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    hashes = {str(p.relative_to(out)).replace("\\", "/"): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in sorted(out.rglob("*")) if p.is_file()}
    write_json(out / "manifest.json", {"script_sha256": script_hash, "files_sha256": hashes})
    print(json.dumps(summary, ensure_ascii=False, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, help="New directory for deterministic scientific outputs and runtime log")
    parser.add_argument("--config", type=Path, help="JSON partial overrides of explicit synthetic defaults")
    parser.add_argument("--self-test", action="store_true", help="Run numerical invariants; combine with --out to run the full study")
    parser.add_argument("--write-example", type=Path, help="Write complete parameter JSON")
    args = parser.parse_args()
    cfg = dict(DEFAULT)
    if args.config:
        overrides = json.loads(args.config.read_text(encoding="utf-8"))
        unknown = set(overrides) - set(cfg)
        if unknown:
            parser.error(f"Unknown configuration keys: {unknown}")
        cfg.update(overrides)
    validate_config(cfg)
    if args.write_example:
        if args.write_example.exists():
            raise FileExistsError(args.write_example)
        write_json(args.write_example, cfg)
    if args.self_test:
        print(json.dumps(numeric_checks(cfg), indent=2))
    if args.out:
        run(cfg, args.out)
    elif not (args.self_test or args.write_example):
        parser.error("Provide --out NEW_DIRECTORY, --self-test, or --write-example PATH")


if __name__ == "__main__":
    main()
