#!/usr/bin/env python3
"""Task 4: reproducible seven-state PBPK research model and dose screening.

Python >=3.10; numpy >=1.24, scipy >=1.10, matplotlib >=3.7.
Run from repo root: python projects/task04_pbpk/run_task4_pbpk_pharmacokinetics_dose_prediction.py --out work/task4_reproduction
Export inputs: --write-example compound_example.json
Use inputs: --config projects/task04_pbpk/compound_example.json --out work/task4_custom
Optional: --self-test, --method BDF, --git-sync.

The default inputs are SYNTHETIC. Outputs are research scenarios, not clinical
dose recommendations. Seven drug-amount states include a GI luminal depot,
NOT a perfused gut wall. Four additional states are bookkeeping integrals.
The script generates all reports, tables, input templates and figures itself.
No network, credentials, hidden datasets or automatic package installation.
"""
from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import math
import platform
import subprocess
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.integrate import solve_ivp, trapezoid
from scipy.linalg import expm
from scipy.optimize import brentq, root

NAMES = ("blood", "liver", "gut_depot", "kidney", "brain", "lung", "rest")
BLOOD, LIVER, GUT, KIDNEY, BRAIN, LUNG, REST = range(7)
HEP, REN, FEC, AUC = range(7, 11)
SYSTEMIC = ("liver", "kidney", "brain", "rest")
TISSUES = ("liver", "kidney", "brain", "lung", "rest")
IDX = {n: i for i, n in enumerate(NAMES)}
COMMIT_MESSAGE = "feat(task4): whole-body PBPK pharmacokinetic model, IVIVE scaling & clinical dose titration"


@dataclass
class Config:
    compound_name: str = "SYNTHETIC_LEAD_DEMO"
    synthetic_inputs: bool = True
    mw_g_mol: float = 450.0
    clogp: float = 3.0
    ionization: str = "base"  # neutral, acid, base, ampholyte
    pka_acid: float | None = None
    pka_base: float | None = 7.8
    fu_p: float = 0.10
    rb: float = 1.0
    clint_mic_ul_min_mg: float = 8.0
    fu_mic: float = 1.0
    mppgl_mg_g: float = 45.0
    liver_mass_g: float = 1800.0
    ka_h: float = 1.0
    fa: float = 0.90
    fg: float = 1.0
    gfr_l_h: float = 7.5
    km_unbound_mg_l: float | None = None
    ic90_nm: float = 50.0
    toxic_total_plasma_mg_l: float = 5.0
    coverage_goal_pct: float = 90.0
    dose_min_mg: float = 10.0
    dose_max_mg: float = 800.0
    dose_step_mg: float = 1.0
    target_tissue: str = "rest"
    blood_volume_l: float = 5.0
    volumes_l: dict = field(default_factory=lambda: {
        "liver": 1.8, "kidney": 0.31, "brain": 1.4,
        "lung": 0.50, "rest": 60.99})
    flows_l_h: dict = field(default_factory=lambda: {
        "liver": 90.0, "kidney": 66.0, "brain": 45.0, "rest": 159.0})
    # Effective volume fractions (neutral lipid, phospholipid, water).
    # Explicit scenario assumptions, not a validated human composition database.
    tissue_fractions: dict = field(default_factory=lambda: {
        "liver": [0.0138, 0.0303, 0.705],
        "kidney": [0.0121, 0.0240, 0.783],
        "brain": [0.0391, 0.0533, 0.770],
        "lung": [0.0030, 0.0090, 0.811],
        "rest": [0.0500, 0.0150, 0.650]})
    plasma_fractions: list = field(default_factory=lambda: [0.00147, 0.00083, 0.96])
    interstitial_protein_ratio: float = 0.5
    kp_overrides: dict = field(default_factory=dict)
    rtol: float = 1e-8
    atol: float = 1e-10
    steady_rtol: float = 1e-7
    max_steady_cycles: int = 2000
    max_tail_h: float = 100000.0

    def validate(self):
        numeric = ("mw_g_mol", "fu_p", "rb", "fu_mic", "mppgl_mg_g",
                   "liver_mass_g", "ka_h", "ic90_nm", "toxic_total_plasma_mg_l",
                   "dose_min_mg", "dose_max_mg", "dose_step_mg", "blood_volume_l",
                   "rtol", "atol", "steady_rtol", "max_tail_h")
        for name in numeric:
            v = getattr(self, name)
            if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v) or v <= 0:
                raise ValueError(f"{name} must be a finite positive number")
        for name in ("clint_mic_ul_min_mg", "gfr_l_h", "interstitial_protein_ratio"):
            v = getattr(self, name)
            if not isinstance(v, (int, float)) or not math.isfinite(v) or v < 0:
                raise ValueError(f"{name} must be finite and nonnegative")
        for name in ("fa", "fg"):
            v = getattr(self, name)
            if not isinstance(v, (int, float)) or not math.isfinite(v) or not 0 <= v <= 1:
                raise ValueError(f"{name} must be in [0, 1]")
        if not 0.01 <= self.fu_p <= 0.20 or not 0.6 <= self.rb <= 1.4 or self.fu_mic > 1:
            raise ValueError("Required ranges: fu_p 0.01..0.20, Rb 0.6..1.4, fu_mic (0,1]")
        if not math.isfinite(self.clogp) or not -5 <= self.clogp <= 10:
            raise ValueError("clogp must be finite and between -5 and 10")
        if self.ionization not in {"neutral", "acid", "base", "ampholyte"}:
            raise ValueError("ionization must be neutral, acid, base or ampholyte")
        for name in ("pka_acid", "pka_base"):
            v = getattr(self, name)
            if v is not None and (not math.isfinite(v) or not -5 <= v <= 20):
                raise ValueError(f"{name} must be null or finite in [-5,20]")
        if self.ionization in {"acid", "ampholyte"} and self.pka_acid is None:
            raise ValueError("pka_acid is required for acid/ampholyte")
        if self.ionization in {"base", "ampholyte"} and self.pka_base is None:
            raise ValueError("pka_base is required for base/ampholyte")
        if self.km_unbound_mg_l is not None:
            if not math.isfinite(self.km_unbound_mg_l) or self.km_unbound_mg_l <= 0:
                raise ValueError("km_unbound_mg_l must be null or finite and positive")
            if self.clint_mic_ul_min_mg <= 0:
                raise ValueError("Saturable hepatic metabolism requires positive CLint")
        if not self.clint_mic_ul_min_mg and not self.gfr_l_h:
            raise ValueError("At least one clearance pathway is required for finite AUC and steady state")
        if not isinstance(self.coverage_goal_pct, (int, float)) or not 0 < self.coverage_goal_pct <= 100:
            raise ValueError("coverage_goal_pct must be in (0,100]")
        if not 10 <= self.dose_min_mg <= self.dose_max_mg <= 800:
            raise ValueError("Per-administration dose bounds must satisfy 10 <= min <= max <= 800 mg")
        if self.target_tissue not in TISSUES:
            raise ValueError("target_tissue must be a modeled tissue")
        if type(self.synthetic_inputs) is not bool or not isinstance(self.compound_name, str):
            raise ValueError("synthetic_inputs must be boolean and compound_name must be text")
        if type(self.max_steady_cycles) is not int or self.max_steady_cycles < 2:
            raise ValueError("max_steady_cycles must be an integer >= 2")
        for mapping, keys in ((self.volumes_l, TISSUES), (self.flows_l_h, SYSTEMIC)):
            if set(mapping) != set(keys):
                raise ValueError(f"Physiology keys must be exactly {keys}")
            if any(not isinstance(v, (int, float)) or not math.isfinite(v) or v <= 0 for v in mapping.values()):
                raise ValueError("Organ volumes and blood flows must be positive and finite")
        if set(self.tissue_fractions) != set(TISSUES):
            raise ValueError("tissue_fractions must specify each tissue")
        for fractions in [self.plasma_fractions, *self.tissue_fractions.values()]:
            if len(fractions) != 3 or any(not math.isfinite(x) or x < 0 for x in fractions) or not 0 < sum(fractions) <= 1:
                raise ValueError("Fractions must contain nonnegative lipid, phospholipid, water with sum (0,1]")
        for k, v in self.kp_overrides.items():
            if k not in TISSUES or not math.isfinite(v) or v <= 0:
                raise ValueError("Kp overrides require recognized tissues and finite positive values")


def partition_coefficients(c: Config):
    """Original Poulin-Theil composition form with explicit logD approximation.

    Kp = [(D*(nl_t+.3*ph_t)+w_t+.7*ph_t) /
          (D*(nl_p+.3*ph_p)+w_p+.7*ph_p)] * fu_p/fu_t.
    fu_t = 1/[1+protein_ratio*(1/fu_p-1)].
    D7.4 uses a monoprotic Henderson-Hasselbalch neutral-fraction estimate.
    This ionization adaptation is NOT the full Rodgers-Rowland algorithm and
    is unsuitable for confident predictions of strong bases/BBB transport.
    Ampholyte mode is a simple two-ionizable-group approximation, not a
    microstate model. Measured Kp values can override each prediction.
    """
    ion = 1.0
    if c.ionization in {"acid", "ampholyte"}:
        ion += 10 ** (7.4 - c.pka_acid)
    if c.ionization in {"base", "ampholyte"}:
        ion += 10 ** (c.pka_base - 7.4)
    d = 10 ** c.clogp / ion
    fut = 1 / (1 + c.interstitial_protein_ratio * (1 / c.fu_p - 1))

    def affinity(f):
        nl, ph, w = f
        return d * (nl + 0.3 * ph) + w + 0.7 * ph

    denominator = affinity(c.plasma_fractions)
    if denominator <= 0:
        raise ValueError("Plasma composition has zero partition capacity")
    kp = {n: affinity(c.tissue_fractions[n]) / denominator * c.fu_p / fut for n in TISSUES}
    kp.update(c.kp_overrides)
    if any(not math.isfinite(x) or x <= 0 for x in kp.values()):
        raise ValueError("Kp estimates must be finite and positive; review composition")
    return kp, {"method": "Poulin-Theil composition form with logD7.4 approximation",
                "logd74": math.log10(d), "estimated_fu_t": fut,
                "neutral_fraction74": 1 / ion}


class PBPK:
    def __init__(self, config: Config, method="Radau"):
        config.validate()
        self.c, self.method = config, method
        self.kp, self.partition_info = partition_coefficients(config)
        self.q = sum(config.flows_l_h.values())
        # uL/min/mg * mg/g * g * 60 min/h / 1e6 uL/L.
        self.clint = (config.clint_mic_ul_min_mg / config.fu_mic *
                      config.mppgl_mg_g * config.liver_mass_g * 60e-6)
        self.clrenal = config.gfr_l_h * config.fu_p  # plasma-referenced clearance
        qh = config.flows_l_h["liver"]
        x = config.fu_p / config.rb * self.clint
        self.clh_blood = qh * x / (qh + x)
        self.fh = 1 - self.clh_blood / qh
        self.linear = config.km_unbound_mg_l is None
        # Linear limit, also used for eigenvalues and nonlinear low-dose tails.
        self.a = np.column_stack([self.derivative(np.eye(7)[i], force_linear=True)[:7] for i in range(7)])
        self.cp_factor = 1 / (config.blood_volume_l * config.rb)
        if np.max(np.linalg.eigvals(self.a).real) >= 0:
            raise ValueError("Model is not asymptotically stable; review clearance and physiology")

    def derivative(self, y, force_linear=False):
        c = self.c
        m = y[:7]
        cp = m[BLOOD] / (c.blood_volume_l * c.rb)
        cb = c.rb * cp
        ven = {n: c.rb * m[IDX[n]] / (c.volumes_l[n] * self.kp[n]) for n in TISSUES}
        arterial = ven["lung"]
        dm = np.zeros(11)
        # Lung is in SERIES with the systemic circuit, never parallel at Qc.
        dm[BLOOD] = sum(c.flows_l_h[n] * ven[n] for n in SYSTEMIC) - self.q * cb
        dm[LUNG] = self.q * (cb - arterial)
        for n in SYSTEMIC:
            dm[IDX[n]] = c.flows_l_h[n] * (arterial - ven[n])
        absorbed = c.ka_h * m[GUT]
        dm[GUT] = -absorbed
        dm[LIVER] += c.fa * c.fg * absorbed
        cu_liver = c.fu_p * ven["liver"] / c.rb
        if self.linear or force_linear:
            hepatic = self.clint * cu_liver
        else:
            km = c.km_unbound_mg_l
            # Vmax = CLint*Km; CLint denotes the low-concentration slope.
            hepatic = self.clint * km * max(cu_liver, 0.0) / (km + max(cu_liver, 0.0))
        renal = self.clrenal * cp
        dm[LIVER] -= hepatic
        dm[BLOOD] -= renal  # lumped renal sink; no renal clearance double-counting
        dm[HEP], dm[REN], dm[FEC], dm[AUC] = hepatic, renal, (1 - c.fa * c.fg) * absorbed, cp
        return dm

    def segment(self, state, duration, rtol=None, atol=None):
        y0 = np.zeros(11)
        y0[:len(state)] = state
        sol = solve_ivp(lambda t, y: self.derivative(y), (0, duration), y0,
                        method=self.method, rtol=rtol or self.c.rtol,
                        atol=atol or self.c.atol, dense_output=True)
        if not sol.success or not np.all(np.isfinite(sol.y)):
            raise RuntimeError(f"ODE solver failed: {sol.message}")
        if np.min(sol.y[:7]) < -max(1e-7, np.max(y0[:7]) * 1e-8):
            raise RuntimeError("Materially negative drug amounts detected")
        conserved = sol.y[:10].sum(axis=0)
        error = np.max(np.abs(conserved - y0[:10].sum()))
        if error > 2e-6 * max(1.0, abs(y0[:10].sum())):
            raise RuntimeError(f"Mass balance failed: {error:g} mg")
        return sol

    def cp(self, state):
        return state[BLOOD] * self.cp_factor

    def free_nm(self, state):
        return self.cp(state) * self.c.fu_p * 1e6 / self.c.mw_g_mol

    def dose_state(self, dose, route):
        state = np.zeros(7)
        state[BLOOD if route == "iv" else GUT] = dose
        return state

    def ivive(self):
        return {"clint_unbound_l_h": self.clint, "fu_b": self.c.fu_p / self.c.rb,
                "clh_blood_l_h": self.clh_blood, "clh_plasma_l_h": self.clh_blood * self.c.rb,
                "eh_low_concentration": 1 - self.fh, "fh_low_concentration": self.fh,
                "fa_fg_fh_pct_linear": 100 * self.c.fa * self.c.fg * self.fh,
                "clrenal_plasma_l_h": self.clrenal, "cardiac_output_l_h": self.q,
                "vmax_mg_h": None if self.linear else self.clint * self.c.km_unbound_mg_l}


def extrema(sol, model, duration):
    """Bracket stationary points on a mixed early/log/linear mesh, then refine."""
    t = np.unique(np.r_[0, duration, np.geomspace(1e-8, duration, 900), np.linspace(0, duration, 1201)])
    def slope(x):
        return model.derivative(sol.sol(x))[BLOOD] * model.cp_factor
    slopes = np.array([slope(x) for x in t])
    stationary = [0., float(duration)]
    for i in np.flatnonzero(slopes[:-1] * slopes[1:] < 0):
        stationary.append(brentq(slope, t[i], t[i+1], xtol=1e-11))
    candidates = np.unique(stationary)
    values = model.cp(sol.sol(candidates))
    imax, imin = int(np.argmax(values)), int(np.argmin(values))
    return {"cmax_mg_l": float(values[imax]), "tmax_h": float(candidates[imax]),
            "cmin_mg_l": float(values[imin]), "tmin_h": float(candidates[imin]),
            "trough_predose_mg_l": float(model.cp(sol.sol(duration)))}, candidates


def coverage(sol, model, tau, dose_scale=1.0, turns=None):
    """Integrate time above threshold using refined crossings, not point counts.

    At periodic steady state, a 12h interval repeated twice has exactly the
    same covered fraction as the full 24h day. Include all stationary points
    to catch short excursions near maxima/minima.
    """
    if turns is None:
        _, turns = extrema(sol, model, tau)
    threshold = model.c.ic90_nm
    def f(t):
        return float(model.free_nm(sol.sol(t))) * dose_scale - threshold
    points = np.unique(np.r_[turns, np.linspace(0, tau, 257)])
    values = np.array([f(t) for t in points])
    crossings = [0., float(tau)]
    for i in np.flatnonzero(values[:-1] * values[1:] < 0):
        crossings.append(brentq(f, points[i], points[i+1], xtol=1e-10))
    crossings.extend(points[values == 0].tolist())
    crossings = sorted(set(crossings))
    covered = sum(b - a for a, b in zip(crossings[:-1], crossings[1:]) if f((a+b)/2) > 0)
    return float(100 * covered / tau)


def single_dose(model, route, dose=100.):
    initial = model.dose_state(dose, route)
    sol = model.segment(initial, 48.)
    metrics, _ = extrema(sol, model, 48.)
    # AUC-infinity is exact for the linear model, including its entire tail.
    state48 = sol.y[:7, -1]
    if model.linear:
        total_auc = float(model.cp(np.linalg.solve(-model.a, initial)))
        tail_auc = float(model.cp(np.linalg.solve(-model.a, state48)))
        tail_residual_fraction = 0.0
        tail_end = None
    else:
        state, elapsed = sol.y[:, -1].copy(), 48.
        while state[:7].sum() > dose * 1e-9:
            if elapsed >= model.c.max_tail_h:
                raise RuntimeError("Nonlinear AUC tail did not converge before max_tail_h")
            dt = min(max(48., elapsed), model.c.max_tail_h - elapsed)
            tail = model.segment(state, dt)
            state, elapsed = tail.y[:, -1], elapsed + dt
        # Remaining low-concentration tail is negligible; retain an explicit
        # linear-limit correction and expose its fraction in the outputs.
        correction = float(model.cp(np.linalg.solve(-model.a, state[:7])))
        total_auc = float(state[AUC] + correction)
        tail_auc = total_auc - float(sol.y[AUC, -1])
        tail_residual_fraction = float(state[:7].sum() / dose)
        tail_end = elapsed
    # Disposition half-life excludes oral absorption; oral apparent terminal
    # half-life includes flip-flop absorption. Both are asymptotic model rates.
    physical = [i for i in range(7) if i != GUT]
    disposition_rate = -float(np.max(np.linalg.eigvals(model.a[np.ix_(physical, physical)]).real))
    apparent_rate = min(disposition_rate, model.c.ka_h) if route == "oral" else disposition_rate
    mass_error = float(np.max(abs(sol.y[:10].sum(axis=0) - dose)))
    metrics.update({"dose_mg": dose, "route": route, "auc_0_48_mg_h_l": float(sol.y[AUC, -1]),
                    "auc_0_inf_mg_h_l": total_auc, "auc_tail_after_48_pct": 100 * tail_auc / total_auc if total_auc else None,
                    "disposition_terminal_half_life_h": math.log(2) / disposition_rate,
                    "apparent_terminal_half_life_h": math.log(2) / apparent_rate,
                    "half_life_definition": "asymptotic low-concentration model eigenmode; not a 48h regression",
                    "nonlinear_tail_end_h": tail_end, "nonlinear_tail_residual_dose_fraction": tail_residual_fraction,
                    "mass_balance_max_error_mg": mass_error})
    # Reconcile tail AUC with the separately integrated AUC for linear mode.
    if model.linear and not math.isclose(metrics["auc_0_48_mg_h_l"] + tail_auc, total_auc, rel_tol=2e-6, abs_tol=1e-7):
        raise RuntimeError("AUC integral and analytic tail disagree")
    return sol, metrics


def steady_state(model, dose, tau):
    jump = model.dose_state(dose, "oral")
    if model.linear or model.c.fa * model.c.fg == 0:
        e = expm(model.a * tau)
        post = np.linalg.solve(np.eye(7) - e, jump)
        if np.min(post) < -1e-7 or not np.all(np.isfinite(post)):
            raise RuntimeError("Invalid periodic steady state")
        sol = model.segment(post, tau)
        residual = float(np.linalg.norm(sol.y[:7, -1] + jump - post, ord=np.inf) / max(1., np.max(post)))
        cycles = 0
    else:
        # With no renal clearance, mean input >= Vmax cannot have steady state.
        if model.clrenal == 0 and dose * model.c.fa * model.c.fg / tau >= model.clint * model.c.km_unbound_mg_l:
            return None, {"status": "no_finite_steady_state", "dose_mg": dose, "tau_h": tau}, None
        # Solve the nonlinear periodic boundary by positive log-space shooting.
        # GI post-dose mass is known exactly and is excluded from root finding.
        # This avoids thousands of transient cycles for long terminal tails.
        e = expm(model.a * tau)
        seed = np.linalg.solve(np.eye(7)-e,jump)
        physical = [i for i in range(7) if i != GUT]
        post = np.maximum(seed,1e-15)
        # Mean-flux approximation provides a close initial guess even when
        # saturation raises exposure far above the low-concentration limit.
        # It is ONLY an initializer; the full time-varying ODE determines SS.
        mean_input=dose*model.c.fa*model.c.fg/tau
        liver_factor=1+model.clrenal/(model.c.flows_l_h["liver"]*model.c.rb)
        km=model.c.km_unbound_mg_l
        vmax=model.clint*km
        if model.clrenal>0:
            def mean_balance(cp):
                cu=model.c.fu_p*cp*liver_factor
                return model.clrenal*cp+vmax*cu/(km+cu)-mean_input
            cp_guess=brentq(mean_balance,0,max(mean_input/model.clrenal,1e-10))
        else:
            cp_guess=km*mean_input/(vmax-mean_input)/model.c.fu_p
        post[BLOOD]=cp_guess*model.c.rb*model.c.blood_volume_l
        for name in TISSUES:
            post[IDX[name]]=cp_guess*model.kp[name]*model.c.volumes_l[name]
        post[LIVER]*=liver_factor
        post[GUT] = dose / (-math.expm1(-model.c.ka_h*tau))
        evaluations = 0
        def shoot(log_mass):
            nonlocal evaluations
            evaluations += 1
            trial=post.copy()
            trial[physical]=np.exp(np.clip(log_mass,-40,40))
            segment=model.segment(trial,tau)
            end=segment.y[:7,-1]
            return np.log(np.maximum(end[physical],1e-30)/trial[physical])
        fit=root(shoot,np.log(post[physical]),method="hybr",
                 options={"xtol":1e-8,"eps":1e-5,"maxfev":180})
        post[physical]=np.exp(np.clip(fit.x,-40,40))
        sol=model.segment(post,tau)
        residual=float(np.linalg.norm(sol.y[:7,-1]+jump-post,ord=np.inf)/max(1.,np.max(post)))
        cycles=0
        if residual > model.c.steady_rtol:
            # Robust fixed-point fallback. A successful solver flag alone is
            # never accepted as evidence of a converged periodic state.
            post=jump.copy()
            for cycles in range(1,model.c.max_steady_cycles+1):
                sol=model.segment(post,tau)
                next_post=sol.y[:7,-1]+jump
                residual=float(np.linalg.norm(next_post-post,ord=np.inf)/max(1.,np.max(next_post)))
                if residual <= model.c.steady_rtol:
                    break
                post=next_post
            else:
                raise RuntimeError(f"Steady state did not converge for {dose:g} mg q{tau:g}h")
    if residual > max(model.c.steady_rtol, 3e-7):
        raise RuntimeError(f"Periodic boundary condition failed: {residual:g}")
    metrics, turns = extrema(sol, model, tau)
    metrics.update({"status": "converged", "dose_mg": dose, "tau_h": tau,
                    "daily_dose_mg": dose * 24 / tau, "periodic_residual_relative": residual,
                    "iterations": cycles, "auc_tau_mg_h_l": float(sol.y[AUC, -1]),
                    "shooting_evaluations": evaluations if not model.linear and model.c.fa*model.c.fg>0 else 0,
                    "coverage_24h_pct": coverage(sol, model, tau, turns=turns)})
    return sol, metrics, turns


def multidose(model, dose=100., tau=12., end=168.):
    times, states, intervals = [], [], []
    state = np.zeros(11)
    # [0,end): 14 doses at 0..156h; 168h is PRE-dose.
    for start in np.arange(0., end, tau):
        state[GUT] += dose
        duration = min(tau, end - start)
        sol = model.segment(state, duration)
        t = np.linspace(0., duration, max(101, int(duration*30)+1))
        times.append(t + start)
        states.append(sol.sol(t))  # duplicate boundary times preserve depot jumps
        intervals.append(sol)
        state = sol.y[:, -1].copy()
    y, t = np.hstack(states), np.concatenate(times)
    count = len(intervals)
    if not math.isclose(float(state[:10].sum()), count*dose, rel_tol=2e-6):
        raise RuntimeError("Multi-dose administered mass is not conserved")
    return t, y, intervals


def optimize_doses(model):
    c = model.c
    rows, choices, curves = [], [], {}
    for tau, label in ((24., "QD"), (12., "BID")):
        if model.linear:
            unit_sol, unit_metrics, turns = steady_state(model, 1., tau)
        cache = {}
        def evaluate(dose):
            dose = float(dose)
            if dose in cache:
                return cache[dose]
            if model.linear:
                cov = coverage(unit_sol, model, tau, dose_scale=dose, turns=turns)
                peak = unit_metrics["cmax_mg_l"] * dose
                trough = unit_metrics["trough_predose_mg_l"] * dose
                status = "converged"
            else:
                _, m, _ = steady_state(model, dose, tau)
                status = m["status"]
                cov, peak, trough = ((m["coverage_24h_pct"], m["cmax_mg_l"], m["trough_predose_mg_l"])
                                     if status == "converged" else (None, None, None))
            row = {"regimen": label, "dose_mg": dose, "daily_dose_mg": dose*24/tau,
                   "coverage_24h_pct": cov, "cmax_ss_mg_l": peak, "ctrough_ss_mg_l": trough,
                   "feasible": status == "converged" and cov >= c.coverage_goal_pct - 1e-8 and peak < c.toxic_total_plasma_mg_l,
                   "status": status}
            cache[dose] = row
            return row
        # Grid is explicit and includes endpoints. Search integers on a dose
        # lattice by monotone bisection, then verify the immediate predecessor.
        # Monotonicity holds for this positive linear/cooperative saturable
        # elimination model; it must not be assumed for arbitrary extensions.
        n = int(math.floor((c.dose_max_mg - c.dose_min_mg) / c.dose_step_mg + 1e-9))
        if n + 1 > 1_000_000:
            raise ValueError("Dose lattice is too large; increase dose_step_mg")
        lattice = c.dose_min_mg + c.dose_step_mg * np.arange(n+1)
        lo, hi = 0, len(lattice)-1
        while lo < hi:
            mid = (lo+hi)//2
            row = evaluate(lattice[mid])
            # Failure of a finite SS is beyond saturation; search lower doses.
            reached = row["status"] != "converged" or row["coverage_24h_pct"] >= c.coverage_goal_pct - 1e-8
            if reached:
                hi = mid
            else:
                lo = mid + 1
        candidate = evaluate(lattice[lo])
        chosen = copy.deepcopy(candidate) if candidate["feasible"] else None
        if chosen and lo > 0 and evaluate(lattice[lo-1])["feasible"]:
            raise RuntimeError("Dose search minimality check failed")
        if chosen:
            chosen["minimum_definition"] = f"minimum on {c.dose_step_mg:g} mg per-administration lattice starting at {c.dose_min_mg:g} mg"
            choices.append(chosen)
        base_grid = np.linspace(c.dose_min_mg, c.dose_max_mg, 100 if model.linear else 24)
        local_grid = np.clip(lattice[lo] + c.dose_step_mg*np.arange(-2,3), c.dose_min_mg,c.dose_max_mg)
        fine_grid = np.geomspace(c.dose_min_mg,c.dose_max_mg,100) if model.linear else []
        grid = np.unique(np.r_[base_grid, fine_grid, local_grid])
        curve = [evaluate(d) for d in grid]
        curves[label] = curve
        rows.extend(cache[d] for d in sorted(cache))
    selected = min(choices, key=lambda r: (r["daily_dose_mg"], r["cmax_ss_mg_l"])) if choices else None
    return {"selected_research_regimen": selected, "minimum_by_regimen": choices,
            "objective": "minimum total daily dose; ties resolved by lower peak",
            "dose_bounds_are_per_administration": True,
            "toxicity_basis": "total mixed-venous plasma concentration (mg/L)",
            "coverage_basis": "free mixed-venous plasma (nM), periodic steady state",
            "clinical_recommendation": None}, rows, curves


def write_csv(path, rows):
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path, obj):
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8", newline="\n")


def profile_rows(t, y, model):
    result = []
    for j, time in enumerate(t):
        row = {"time_h": float(time), "plasma_mg_l": float(model.cp(y[:, j])),
               "free_plasma_nm": float(model.free_nm(y[:, j]))}
        row.update({f"{n}_mg_l": float(y[IDX[n], j] / model.c.volumes_l[n]) for n in TISSUES})
        row.update({f"amount_{n}_mg": float(y[i,j]) for i, n in enumerate(NAMES)})
        row.update({"hepatic_eliminated_mg": float(y[HEP,j]), "renal_eliminated_mg": float(y[REN,j]),
                    "presystemic_loss_mg": float(y[FEC,j]), "auc_to_time_mg_h_l": float(y[AUC,j])})
        result.append(row)
    return result


def make_figures(out, model, iv, oral, metrics, mt, my, ss, curves):
    folder = out / "figures_task4"
    folder.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.titleweight": "bold", "savefig.dpi": 300})
    stamp = "SYNTHETIC INPUTS | RESEARCH SCENARIO" if model.c.synthetic_inputs else "UNVALIDATED RESEARCH SCENARIO"
    def save(fig, name):
        fig.text(0.5, 0.012, stamp, ha="center", fontsize=8, color="#7b4855")
        fig.tight_layout(rect=(0, .045, 1, .95))
        fig.savefig(folder/name, dpi=300, facecolor="white")
        plt.close(fig)
    t = np.unique(np.r_[0, np.geomspace(1e-5, 48, 1800), np.linspace(0,48,1201)])
    fig, axes = plt.subplots(1,2, figsize=(12,4.9))
    for route, sol, color in (("iv", iv, "#2458a6"), ("oral", oral, "#d96629")):
        m = metrics[route]
        label = (f"{route.upper()}: Cmax {m['cmax_mg_l']:.3g} mg/L; Tmax {m['tmax_h']:.3g} h\n"
                 f"AUC0-inf {m['auc_0_inf_mg_h_l']:.3g} mg h/L")
        cp = model.cp(sol.sol(t))
        for ax in axes:
            ax.plot(t, np.where(cp>0, cp, np.nan), color=color, lw=2, label=label)
            ax.scatter([m["tmax_h"]], [m["cmax_mg_l"]], color=color, s=25, zorder=4)
    axes[1].set_yscale("log")
    visible = np.r_[model.cp(iv.sol(t[t>=.05])), model.cp(oral.sol(t[t>=.05]))]
    positive = visible[visible>0]
    if positive.size:
        axes[1].set_ylim(max(1e-12,positive.min()/3), metrics["iv"]["cmax_mg_l"]*1.6)
    inset=axes[0].inset_axes([.26,.24,.68,.42])
    for sol,color in ((iv,"#2458a6"),(oral,"#d96629")):
        inset.plot(t[t<=12],model.cp(sol.sol(t[t<=12])),color=color,lw=1.3)
    inset.set(xlim=(0,12),ylim=(0,max(metrics["oral"]["cmax_mg_l"]*1.25,.01)),
              title="0-12 h: expanded concentration scale")
    inset.tick_params(labelsize=8); inset.title.set_fontsize(8); inset.grid(alpha=.15)
    for ax, title in zip(axes, ("Linear scale", "Semi-log scale")):
        ax.set(xlabel="Time (h)", ylabel="Total plasma concentration (mg/L)", title=title, xlim=(0,48))
        ax.grid(alpha=.18); ax.legend(fontsize=8, loc="upper right")
    fig.suptitle("Single 100 mg dose: IV versus oral", fontsize=15)
    save(fig, "fig1_pbpk_plasma_iv_vs_oral.png")

    fig, ax = plt.subplots(figsize=(10,5.5))
    t = np.linspace(0,48,1801); y = oral.sol(t)
    ax.plot(t, model.cp(y), label="Plasma", lw=2.5, color="#253858")
    organs = list(dict.fromkeys(["liver", "brain", "kidney", model.c.target_tissue]))
    for n in organs:
        label = n.title() + (" (target proxy)" if n == model.c.target_tissue else "")
        ax.plot(t, y[IDX[n]]/model.c.volumes_l[n], label=label, lw=1.8)
    ax.set(xlabel="Time (h)", ylabel="Total concentration (mg/L)", xlim=(0,48), title="Tissue distribution after oral 100 mg")
    ax.legend(); ax.grid(alpha=.18)
    save(fig, "fig2_tissue_distribution_biodistribution.png")

    fig, axes = plt.subplots(2,1,figsize=(11,7), gridspec_kw={"height_ratios":[1.4,1]})
    axes[0].plot(mt,model.cp(my),color="#2458a6",lw=1.6,label="Day 1-7 profile")
    for dose_time in np.arange(0,168,12):
        axes[0].axvline(dose_time,alpha=.09,color="black",lw=.7)
    if ss[1]["status"] == "converged":
        sm = ss[1]
        axes[0].axhline(sm["cmax_mg_l"],color="#d96629",ls="--",label=f"True SS peak {sm['cmax_mg_l']:.3g} mg/L")
        axes[0].axhline(sm["trough_predose_mg_l"],color="#278568",ls="--",label=f"True SS pre-dose trough {sm['trough_predose_mg_l']:.3g} mg/L")
        st = np.linspace(0,12,601)
        axes[1].plot(st,model.free_nm(ss[0].sol(st)),lw=2,label="True steady-state interval")
    axes[1].plot(mt[mt>=156]-156, model.free_nm(my[:,mt>=156]),lw=1.8,ls=":",label="Last interval on day 7")
    axes[1].axhline(model.c.ic90_nm, color="#b44259",ls="--",label=f"IC90 = {model.c.ic90_nm:g} nM")
    axes[0].set(xlabel="Time (h)",ylabel="Total plasma (mg/L)",xlim=(0,168),title="100 mg BID: 14 doses at 0-156 h; endpoint 168 h is pre-dose")
    axes[1].set(xlabel="Time since last dose (h)",ylabel="Free plasma (nM)",xlim=(0,12))
    for ax in axes: ax.legend(fontsize=9); ax.grid(alpha=.18)
    save(fig,"fig3_multidose_steady_state_regimen.png")

    fig, axes = plt.subplots(1,2,figsize=(12,5.4))
    dose_inset=axes[0].inset_axes([.15,.30,.46,.36])
    for label,color in (("QD","#2458a6"),("BID","#d96629")):
        curve=curves[label]
        x=np.array([r["daily_dose_mg"] for r in curve])
        cov=np.array([r["coverage_24h_pct"] if r["coverage_24h_pct"] is not None else np.nan for r in curve])
        peak=np.array([r["cmax_ss_mg_l"] if r["cmax_ss_mg_l"] is not None else np.nan for r in curve])
        feasible=np.array([r["feasible"] for r in curve])
        axes[0].plot(x,cov,label=label,color=color,lw=2)
        axes[0].fill_between(x,0,100,where=feasible,color=color,alpha=.10,label=f"{label} feasible at plotted doses")
        axes[1].plot(x,peak,label=label,color=color,lw=2)
        dose_inset.plot(x,cov,color=color,lw=1.4)
    dose_inset.axhline(model.c.coverage_goal_pct,color="#278568",ls="--",lw=1)
    dose_inset.set(xlim=(0,min(160,model.c.dose_max_mg*2)),ylim=(0,102),title="Low-dose detail (mg/day)")
    dose_inset.tick_params(labelsize=8); dose_inset.title.set_fontsize(8); dose_inset.grid(alpha=.15)
    axes[0].axhline(model.c.coverage_goal_pct,color="#278568",ls="--",label="Coverage requirement")
    axes[1].axhline(model.c.toxic_total_plasma_mg_l,color="#b44259",ls="--",label="Assumed toxicity threshold")
    axes[0].set(ylabel="Time above free IC90 per 24 h (%)",ylim=(0,102),title="Coverage by dose and regimen")
    axes[1].set(ylabel="Steady-state total plasma Cmax (mg/L)",title="Peak-exposure constraint")
    for ax in axes: ax.set_xlabel("Total daily oral dose (mg/day)"); ax.grid(alpha=.18); ax.legend(fontsize=8)
    fig.suptitle("Dose screening at periodic steady state (10-800 mg per administration)",fontsize=13)
    save(fig,"fig4_dose_titration_target_coverage.png")


def self_tests():
    """Independent invariants, limiting cases, and solver cross-checks."""
    checks = []
    def check(name, condition):
        if not condition:
            raise AssertionError(name)
        checks.append({"test":name,"passed":True})
    c=Config(); m=PBPK(c)
    check("IVIVE units: 8*45*1800*60e-6 = 38.88 L/h", math.isclose(m.clint,38.88))
    rng=np.random.default_rng(401)
    state=rng.uniform(0,100,11)
    d=m.derivative(state)
    check("Analytic mass derivative including losses is zero", abs(d[:10].sum())<1e-9)
    sol=m.segment(m.dose_state(100,"oral"),48)
    exact=expm(m.a*48)@m.dose_state(100,"oral")
    check("Radau agrees with independent matrix exponential", np.allclose(sol.y[:7,-1],exact,rtol=2e-6,atol=1e-8))
    b=PBPK(c,"BDF").segment(m.dose_state(100,"oral"),48)
    check("BDF agrees with Radau", np.allclose(b.y[:7,-1],sol.y[:7,-1],rtol=3e-6,atol=1e-8))
    for rb in (.6,1.,1.4):
        x=copy.deepcopy(c); x.rb=rb; z=PBPK(x)
        ai=z.cp(np.linalg.solve(-z.a,z.dose_state(100,"iv")))
        ao=z.cp(np.linalg.solve(-z.a,z.dose_state(100,"oral")))
        check(f"AUC ratio matches Fa*Fg*FH at Rb={rb}", math.isclose(ao/ai,x.fa*x.fg*z.fh,rel_tol=1e-10))
    x=copy.deepcopy(c); x.clint_mic_ul_min_mg=0; z=PBPK(x)
    check("Zero hepatic clearance gives FH=1", z.fh==1)
    x=copy.deepcopy(c); x.fa=0; z=PBPK(x)
    zero=z.segment(z.dose_state(100,"oral"),48)
    check("Fa=0 gives zero oral systemic exposure", np.max(np.abs(zero.y[[BLOOD,LIVER,KIDNEY,BRAIN,LUNG,REST]]))<1e-10)
    ss,sm,turns=steady_state(m,100,12)
    check("Periodic state residual is below tolerance",sm["periodic_residual_relative"]<1e-7)
    check("Steady-state AUC per interval equals single-dose AUC infinity",math.isclose(sm["auc_tau_mg_h_l"],m.cp(np.linalg.solve(-m.a,m.dose_state(100,"oral"))),rel_tol=2e-6))
    check("Zero dose has zero target coverage",coverage(ss,m,12,dose_scale=0,turns=turns)==0)
    check("Very high scaled exposure has complete coverage",coverage(ss,m,12,dose_scale=1e8,turns=turns)==100)
    q=PBPK(c); q.clint=0; q.clrenal=0
    closed=q.segment(q.dose_state(100,"iv"),48)
    check("Closed system conserves all active drug",np.allclose(closed.y[:7].sum(axis=0),100,atol=1e-6))
    x=copy.deepcopy(c); x.km_unbound_mg_l=1e8; z=PBPK(x)
    n=z.segment(z.dose_state(100,"oral"),48)
    check("Large-Km nonlinear model approaches linear model",np.allclose(n.y[:7,-1],sol.y[:7,-1],rtol=3e-6,atol=1e-8))
    x=copy.deepcopy(c); x.km_unbound_mg_l=.1; z=PBPK(x)
    check("Nonlinear model mass derivative is zero",abs(z.derivative(state)[:10].sum())<1e-9)
    periodic,pm,_=steady_state(z,100,12)
    check("Nonlinear shooting satisfies each compartment periodic boundary",np.allclose(periodic.y[:7,-1]+z.dose_state(100,"oral"),periodic.y[:7,0],rtol=2e-6,atol=1e-7))
    check("Nonlinear periodic elimination equals absorbed input",math.isclose(float(periodic.y[HEP,-1]+periodic.y[REN,-1]),100*x.fa*x.fg,rel_tol=2e-6))
    x=copy.deepcopy(c); x.km_unbound_mg_l=.001; x.gfr_l_h=0; z=PBPK(x)
    _,unstable,_=steady_state(z,100,12)
    check("Saturable capacity failure is identified",unstable["status"]=="no_finite_steady_state")
    for key,value in (("mw_g_mol",0),("fu_p",.5),("rb",float("nan")),("fa",1.5),("km_unbound_mg_l",-1)):
        x=copy.deepcopy(c); setattr(x,key,value)
        try: x.validate()
        except ValueError: invalid=True
        else: invalid=False
        check(f"Invalid {key} is rejected",invalid)
    # Independent time quadrature sanity check for partial coverage.
    unit,_,ut=steady_state(m,1.,12)
    cv=lambda dose: coverage(unit,m,12,dose_scale=dose,turns=ut)
    partial=brentq(lambda dose:cv(dose)-50,1,10000,xtol=1e-7)
    dense=np.linspace(0,12,200001)
    fraction=trapezoid((m.free_nm(unit.sol(dense))*partial>c.ic90_nm).astype(float),dense)/12*100
    check("Root-based coverage matches fine independent quadrature",abs(cv(partial)-fraction)<.002)
    return checks


def reports(out, model, results):
    c=model.c; iv=results["single_dose"]["iv"]; po=results["single_dose"]["oral"]
    ss=results["bid_100mg_steady_state"]; day=results["day7_last_interval"]
    opt=results["dose_optimization"]; pick=opt["selected_research_regimen"]
    selected_en=(f"{pick['dose_mg']:g} mg {pick['regimen']} ({pick['daily_dose_mg']:g} mg/day), "
                 f"coverage {pick['coverage_24h_pct']:.3f}%, Cmax {pick['cmax_ss_mg_l']:.4f} mg/L") if pick else "No feasible finite-steady-state regimen on the specified dose lattice."
    selected_zh=(f"每次 {pick['dose_mg']:g} mg，{pick['regimen']}（{pick['daily_dose_mg']:g} mg/日），"
                 f"覆盖率 {pick['coverage_24h_pct']:.3f}%，峰浓度 {pick['cmax_ss_mg_l']:.4f} mg/L") if pick else "指定剂量网格中无满足条件的有限稳态方案。"
    kp_table="\n".join(f"| {n} | {c.volumes_l[n]:g} | {model.q if n=='lung' else c.flows_l_h[n]:g} | {model.kp[n]:.6g} | {','.join(map(str,c.tissue_fractions[n]))} |" for n in TISSUES)
    endpoints=(f"| IV Cmax / Tmax | {iv['cmax_mg_l']:.6g} mg/L / {iv['tmax_h']:.6g} h |\n"
               f"| Oral Cmax / Tmax | {po['cmax_mg_l']:.6g} mg/L / {po['tmax_h']:.6g} h |\n"
               f"| IV AUC 0-inf | {iv['auc_0_inf_mg_h_l']:.6g} mg h/L |\n"
               f"| Oral AUC 0-inf | {po['auc_0_inf_mg_h_l']:.6g} mg h/L |\n"
               f"| Oral apparent terminal half-life | {po['apparent_terminal_half_life_h']:.6g} h |\n"
               f"| AUC-ratio oral F | {results['oral_bioavailability_auc_ratio_pct']:.6g}% |\n"
               f"| Oral AUC after 48 h | {po['auc_tail_after_48_pct']}% |\n"
               f"| Day 7 pre-dose trough | {day['trough_predose_mg_l']:.6g} mg/L |")
    if ss["status"]=="converged":
        endpoints+=(f"\n| True SS peak / pre-dose trough | {ss['cmax_mg_l']:.6g} / {ss['trough_predose_mg_l']:.6g} mg/L |"
                    f"\n| True SS within-interval minimum | {ss['cmin_mg_l']:.6g} mg/L |"
                    f"\n| True SS 24h coverage | {ss['coverage_24h_pct']:.6g}% |"
                    f"\n| Cmax accumulation ratio | {ss['accumulation_ratio']} |"
                    f"\n| Day 7 vs SS state relative error | {day['state_error_vs_ss_relative']:.6g} |"
                    f"\n| Day 7 reached SS (1% state/peak/trough criterion) | {day['within_1pct_of_ss']} |")
    refs="""1. Poulin & Theil (2002), tissue-composition partitioning: https://doi.org/10.1002/jps.10005
2. Poulin & Theil (2002), generic PBPK models: https://doi.org/10.1002/jps.10128
3. FDA (2018), PBPK report format and content: https://www.fda.gov/regulatory-information/search-fda-guidance-documents/physiologically-based-pharmacokinetic-analyses-format-and-content-guidance-industry
"""
    equations=r"""
Let M_B be mixed-venous blood drug mass; Cp=M_B/(V_B Rb).
For tissue i, C_i=M_i/V_i and C_vi=Rb*C_i/Kp_i (venous blood).
C_a=C_v,lung; Qc=sum_i Q_i, i=liver,kidney,brain,rest.

    dM_g/dt = -ka*M_g
    dM_lung/dt = Qc*(Rb*Cp - C_a)
    dM_i/dt = Qi*(C_a-C_vi)                      [kidney, brain, rest]
    dM_L/dt = QH*(C_a-C_vL) + ka*Fa*Fg*M_g - H
    dM_B/dt = sum_i(Qi*C_vi) - Qc*Rb*Cp - CLrenal,p*Cp
    H_linear = CLint,u * fu,p * C_L/Kp,L
    H_saturable = Vmax*Cu,L/(Km+Cu,L), Cu,L=fu,p*C_L/Kp,L
    Vmax = CLint,u*Km; CLrenal,p = GFR*fu,p

Bookkeeping: dE_H/dt=H; dE_R/dt=CLrenal,p*Cp;
dE_pre/dt=(1-Fa*Fg)*ka*M_g; dAUC/dt=Cp.
M_B+M_g+sum(M_tissues)+E_H+E_R+E_pre = cumulative administered dose.

    CLint,u [L/h] = (CLint,mic / fu,mic) * MPPGL * liver_mass_g * 60/1e6
    fu,b = fu,p/Rb
    CLH,b = QH*fu,b*CLint,u / (QH+fu,b*CLint,u)
    EH = CLH,b/QH; FH = 1-EH; F_linear = Fa*Fg*FH
    Cfree [nM] = Cp [mg/L] * fu,p * 1e6 / MW [g/mol]

Kp approximation (common pH 7.4):
    D = 10^cLogP / (1 + acid_term + base_term)
    acid_term = 10^(7.4-pKa_acid); base_term = 10^(pKa_base-7.4)
    [include only terms appropriate to the selected ionization class]
    fu,t = 1/(1 + protein_ratio*(1/fu,p-1))
    A(f) = D*(f_nl+0.3*f_ph)+f_water+0.7*f_ph
    Kp,t = A(tissue)/A(plasma) * fu,p/fu,t

Linear model dM/dt=A*M:
    AUC_0_inf = e_B^T*(-A)^(-1)*M0/(V_B*Rb)
    M_ss,post = (I-exp(A*tau))^(-1)*dose_vector
"""
    en=f"""# Task 4 — PBPK and exposure-based dose screening

## Executive summary
Compound: **{c.compound_name}**. Synthetic inputs: **{c.synthetic_inputs}**.
This is a reproducible, reduced seven-state research PBPK model, not a qualified clinical prediction platform. Input provenance, assumptions, numerical checks and uncertainty must be reviewed before interpreting exposure. No lead-specific measured clearance, tissue Kp values, toxicology dataset or human PK observations was supplied for the delivered demonstration. The configured toxicity threshold is a scenario assumption, not an established safe exposure limit. Setting `synthetic_inputs=false` changes the label only; it does not qualify the model.

The selected **research scenario** is: {selected_en}. A clinical recommended dose is deliberately left null in the machine-readable results. QD is every 24 h and BID every 12 h. Doses refer to active compound mass; salt correction is not inferred.

## Translational rationale and model scope
PBPK joins compound properties, organ capacities, circulation and clearance to connect an administered dose with time-varying exposure. Comparing free plasma exposure with potency can rank formulations or regimens and identify measurements that most affect a development decision. Here, IC90 is treated as a free-concentration threshold. If an assay reports a nominal concentration, binding and assay-to-tissue translation must be resolved separately. Time above IC90 is an exposure surrogate; it is not a measurement of receptor occupancy or clinical efficacy.

The requested seven compartments are blood/plasma, liver, GI depot, kidney, brain, lung and remaining tissues. Blood is an amount in whole blood with a mixed-venous plasma observation obtained using Rb. GI is a lumen depot; gut wall, portal blood and arterial blood do not have separate states. Lung is in series with the systemic organs. Hepatic flow lumps hepatic arterial and portal supply. The renal sink is placed in the central state as specified by the brief; the kidney tissue state is distribution-only. This is a reduced circulation approximation. Remaining tissues use an effective composition and act as a target proxy, not a separately validated target organ.

## Inputs and physiology
All full inputs and units are in `parameters_used.json`; an editable template is `compound_example.json`. MW={c.mw_g_mol:g} g/mol, cLogP={c.clogp:g}, ionization={c.ionization}, acidic/basic pKa={c.pka_acid}/{c.pka_base}, fu,p={c.fu_p:g}, Rb={c.rb:g}, CLint,mic={c.clint_mic_ul_min_mg:g} uL/min/mg, fu,mic={c.fu_mic:g}, ka={c.ka_h:g} h^-1, Fa={c.fa:g}, Fg={c.fg:g}. MW is used for unit conversion; it does not independently determine Kp in this model. Blood volume={c.blood_volume_l:g} L. Tissue volumes and flows are fixed adult scenario values, not individualized allometric predictions; 1,800 g is liver mass, distinct from the 1.8 L liver distribution volume.

| Tissue | Volume L | Flow L/h | Kp tissue:plasma | Fractions: neutral lipid, phospholipid, water |
|---|---:|---:|---:|---|
{kp_table}

Tissue partitioning uses the Poulin–Theil composition form [1], with a clearly labeled Henderson–Hasselbalch logD7.4 substitution for ionization. Protein partitioning uses the configured interstitial/plasma ratio. This adaptation assumes equal pH, passive equilibration and binding surrogates. It is **not** a full Rodgers–Rowland or Berezhkovskiy implementation. It does not model acidic phospholipid binding, lysosomal trapping, active transport or BBB permeability. Tissue fractions are explicit illustrative assumptions. Brain and pooled-rest predictions particularly require empirical verification; measured Kp overrides are supported. Ampholytes use a simplified ionization approximation.

## IVIVE and equations
Microsomal clearance is interpreted as apparent clearance based on total incubation concentration and corrected by fu,mic. Set fu,mic=1 if the supplied clearance is already unbound. The explicit factor 60/1e6 converts uL/min to L/h. At the current inputs, CLint,u={model.clint:.6g} L/h, CLH,b={model.clh_blood:.6g} L/h and low-concentration FH={model.fh:.6g}. FH is hepatic escape, not overall oral availability. Use it once through hepatic metabolism; multiplying the liver input by FH again would double-count first pass.

```text
{equations}
```

The default ODEs are linear. Optional `km_unbound_mg_l` enables saturable metabolism with Vmax=CLint,u*Km. Km must be supplied or explicitly assumed; it cannot be inferred from microsomal CLint alone. Well-stirred clearance then describes only the low-concentration limit. With nonlinear kinetics the equal-dose oral/IV AUC ratio remains an apparent exposure ratio and is not generally the absolute fraction reaching circulation; it can exceed 100%. No dose-linear scaling is used in nonlinear dose screening.

## Numerical methods and results
Radau or BDF integrates every dosing interval separately, with exact state jumps at administration. The 7-day simulation gives 14 oral doses at 0,12,...,156 h; 168 h is pre-dose. Seven amount states and four integral counters are distinguished. Matrix solutions independently check linear integration, integrate the entire AUC tail and solve the periodic steady-state boundary. Nonlinear AUC integrates until residual mass is below 1e-9 of the dose, then adds a low-concentration tail approximation. Finite limits trigger errors rather than silently accepting unconverged results. Reported half-lives are asymptotic low-concentration eigenmode rates; the oral rate allows flip-flop absorption and is not a fitted half-life from the 48 h graph. Cmax is searched on 0-48 h for the single doses.

| Endpoint | Value |
|---|---|
{endpoints}

The day-7 result is tested against a separate periodic steady state; seven days is not assumed sufficient. The pre-dose trough and actual within-interval minimum are both reported because continued absorption can move the minimum away from the dosing boundary. Accumulation uses the peak over the first 12 h interval, not an arbitrarily long single-dose window. Complete metrics, mass errors and periodic residuals are in `results_summary.json`.

## Dose rationale, uncertainty and reproducibility
For each QD/BID regimen, the script finds the minimum on the configured per-administration lattice ({c.dose_min_mg:g}–{c.dose_max_mg:g} mg, step {c.dose_step_mg:g} mg), evaluates free-plasma coverage >= {c.coverage_goal_pct:g}% and enforces **total-plasma** Cmax < {c.toxic_total_plasma_mg_l:g} mg/L. The daily-dose objective can therefore span 10–800 mg/day for QD and 20–1600 mg/day for BID with default bounds. Root-refined threshold crossings determine continuous covered time, including strict equality handling. Regimen minima are compared by total daily dose, then peak. The plotted shaded region is illustrative on the displayed grid; the tabulated selected dose is separately verified on the search lattice.

`sensitivity_analysis.csv` changes fu,p, microsomal clearance, absorption, and all Kp values one at a time. These scenario changes are not confidence intervals, a population simulation or global sensitivity analysis. Key missing determinants include formulation/dissolution, intestinal metabolism evidence, transporter effects, organ disease, interindividual variability, toxicology and observed human PK. Independent clinical data would be needed to evaluate predictive performance and clinical dose selection. FDA guidance [3] distinguishes intended use and the evidence supporting a PBPK report; passing code tests alone does not establish clinical validity.

Run `python projects/task04_pbpk/run_task4_pbpk_pharmacokinetics_dose_prediction.py --self-test --out work/task4_reproduction`. To replace assumptions, export a template with `--write-example work/task4_compound.json`, edit it, then use `--config work/task4_compound.json`. `verification.json` records numerical checks; `manifest.json` contains file hashes and versions. Four figures are saved at 300 DPI. The self-contained script regenerates this report and all artifacts. Task 1 descriptors can be transferred through the JSON template; they do not replace measured binding, clearance or toxicology inputs.

## References
{refs}
"""
    zh=f"""# 任务四：PBPK 与基于暴露的剂量筛选

## 摘要
化合物：**{c.compound_name}**；合成示例输入：**{c.synthetic_inputs}**。
本交付为可复现的简化七状态 PBPK 科研模型，不是经过临床验证的预测平台。未提供指定先导化合物的实测清除率、组织分配实测值、毒理学数据或人体药代观测，因此示例不能产生可靠的人体推荐剂量。毒性阈值为演示假设，不代表已确立的安全上限。将 synthetic_inputs 改成 false 只改变标签，不构成模型验证。

筛选得到的**科研情景方案**为：{selected_zh}。机器可读结果中的 clinical_recommendation 保持为空。QD 表示每 24 小时一次，BID 表示每 12 小时一次；剂量按活性化合物质量计，不自动进行盐型质量换算。

## 药学工程意义与范围
PBPK 将分子性质、器官容量、血流和清除机制连接起来，用于解释给药剂量如何形成动态体内暴露。把游离暴露与药效阈值比较，可用于提出制剂或给药频率的研究假设，并确定下一步应优先测量的参数。这里把 IC90 当作游离浓度阈值；若体外数值是培养体系名义浓度，必须另行校正蛋白结合和测定条件。超过 IC90 的时间比例仅是暴露替代指标，不能等同于受体占有率或临床疗效。

七个药物量状态为血液/血浆、肝、胃肠吸收库、肾、脑、肺、其余组织。中央状态存储全血药物量，再通过 Rb 换算混合静脉血浆浓度。胃肠为腔内药物库，不是灌流肠壁；肠壁、门静脉与动脉血均未单独建室。肺处于体循环串联位置；肝血流汇总门静脉及肝动脉供应。按照任务中的简化肾清除写法，肾排泄从中央室扣除，肾组织室只承担分布，避免重复清除。其余组织是有效成分混合室，默认用作靶组织代理，并非经过验证的独立靶器官。

## 参数与组织分配
完整输入与单位见 parameters_used.json，可编辑模板见 compound_example.json。MW={c.mw_g_mol:g} g/mol，cLogP={c.clogp:g}，电离类型={c.ionization}，酸/碱 pKa={c.pka_acid}/{c.pka_base}，fu,p={c.fu_p:g}，Rb={c.rb:g}，CLint,mic={c.clint_mic_ul_min_mg:g} μL/min/mg，fu,mic={c.fu_mic:g}，ka={c.ka_h:g} h^-1，Fa={c.fa:g}，Fg={c.fg:g}。MW 用于浓度单位换算，不在该分配公式中独立决定 Kp。中央血容量={c.blood_volume_l:g} L；器官体积、血流和组成属于固定成人情景假设，不代表个体化生理参数。肝质量 1,800 g 与分布体积 1.8 L 分开使用。

| 组织 | 体积 L | 血流 L/h | 组织:血浆 Kp | 中性脂质、磷脂、水体积分数 |
|---|---:|---:|---:|---|
{kp_table}

组织分配采用 Poulin–Theil 组成公式[1]，并明确以 Henderson–Hasselbalch 估算的 logD7.4 处理电离；组织蛋白结合采用可配置的间质/血浆蛋白比例。这是等 pH、被动平衡下的近似，**不是完整 Rodgers–Rowland 或 Berezhkovskiy 实现**。未描述酸性磷脂结合、溶酶体捕获、主动转运或血脑屏障通透性。组织组成数值为明确列出的演示假设，脑和混合剩余组织尤其需要实测验证。支持按组织覆盖 Kp；两性化合物仅使用简化电离近似。

## IVIVE 与质量守恒方程
微粒体清除率按基于孵育总浓度的表观清除率解释，除以 fu,mic 后外推；若输入已经是游离清除率，应设 fu,mic=1。乘以 60/10^6 完成 μL/min 到 L/h 的转换。本情景 CLint,u={model.clint:.6g} L/h，肝血清除率 CLH,b={model.clh_blood:.6g} L/h，低浓度肝逃逸率 FH={model.fh:.6g}。FH 仅是肝逃逸率，不是总体口服利用度；肝内代谢已经实现首过，吸收输入不能再次乘 FH。

以下代码块列出完整公式，M 为 mg，浓度为 mg/L，流量及清除率为 L/h，时间为 h。血液流出浓度为 Rb*C_i/Kp_i；这种换算避免将血流与血浆浓度直接相乘。

```text
{equations}
```

默认模型为线性。设置 km_unbound_mg_l 后启用 Michaelis–Menten 饱和代谢，并令 Vmax=CLint,u*Km。仅凭微粒体 CLint 不能推定 Km，必须另行提供或明确假设。此时良好搅拌模型指标只代表低浓度极限；同剂量口服/静脉 AUC 比值应解释为表观暴露比，不一定等于进入循环的绝对剂量比例，也可能超过 100%。非线性剂量筛选逐剂量计算，不使用线性叠加。

## 数值计算与结果
采用 Radau 或 BDF 分段积分，每次给药以精确状态跳跃加入胃肠药物库。0、12、…、156 小时共给药 14 次，168 小时记录给药前值。七个药物量状态之外，另有肝消除、肾消除、首过前损失和 AUC 四个记账积分，不把它们计入生理房室。

线性模型通过矩阵方法独立核验 ODE、计算完整 AUC 尾部，并直接求周期稳态。非线性 AUC 积分至残余质量低于剂量的 10^-9 后，补入低浓度尾部近似；超出迭代或时间上限会报错。半衰期为模型低浓度渐近特征值对应的半衰期，口服指标包含吸收限速的可能，不把 48 小时图形回归值伪装成终末半衰期。单次给药的 Cmax 在 0–48 小时内搜索。

| 指标 | 结果 |
|---|---|
{endpoints}

第七天结果与独立周期稳态比较，不预设七天足够达到稳态。同时报告给药前谷浓度和给药间隔内实际最低浓度，因为吸收延续可能使两者不同。蓄积比使用真正稳态峰浓度除以首次给药后 0–12 小时峰浓度。完整数值、质量守恒误差和周期边界残差见 results_summary.json。

## 剂量筛选、局限与复现
对 QD、BID 分别在每次 {c.dose_min_mg:g}–{c.dose_max_mg:g} mg、步长 {c.dose_step_mg:g} mg 的网格上求最小剂量，要求每 24 小时游离浓度超过 IC90 的时间不少于 {c.coverage_goal_pct:g}%，且**总血浆**峰浓度严格低于 {c.toxic_total_plasma_mg_l:g} mg/L。默认 QD 日剂量范围为 10–800 mg，BID 为 20–1600 mg。阈值交点使用求根精化后按持续时间计算，避免仅统计采样点。两种频率按总日剂量比较，相同日剂量时选择峰值较低者。图中的阴影根据绘图网格展示可行区域，最终候选剂量另在搜索网格中验证。

sensitivity_analysis.csv 提供 fu,p、微粒体清除率、吸收及全部 Kp 的单因素变化结果，用于显示情景敏感性；它不是置信区间、群体变异预测或全局敏感性分析。模型缺少制剂溶出、肠壁代谢实测依据、转运、器官疾病、个体差异、毒理和人体观测。需要外部实验与人体数据检验预测表现，才能支持临床剂量决策。FDA 报告指南[3]强调清楚说明模型用途及支持证据；数值测试通过不等于临床验证。

运行 `python projects/task04_pbpk/run_task4_pbpk_pharmacokinetics_dose_prediction.py --self-test --out work/task4_reproduction` 可复现交付。使用 `--write-example work/task4_compound.json` 导出模板，编辑后通过 `--config work/task4_compound.json` 读取。verification.json 保存检查结果，manifest.json 保存版本和文件校验值，四张图均为 300 DPI。脚本自包含，可重新生成所有表格和中英文报告；Task 1 描述符可通过 JSON 模板导入，但不能替代蛋白结合、清除率或毒理学实测参数。

## 参考资料
{refs}
"""
    (out/"PBPK_DOSE_PREDICTION_REPORT_EN.md").write_text(en.rstrip()+"\n",encoding="utf-8", newline="\n")
    (out/"PBPK_DOSE_PREDICTION_REPORT_ZH.md").write_text(zh.rstrip()+"\n",encoding="utf-8", newline="\n")


def sensitivity(model):
    rows=[]
    scenarios=[("baseline",None,None), ("fu_p_half","fu_p",max(.01,model.c.fu_p*.5)),
               ("fu_p_double","fu_p",min(.2,model.c.fu_p*2)),
               ("clint_half","clint_mic_ul_min_mg",model.c.clint_mic_ul_min_mg*.5),
               ("clint_double","clint_mic_ul_min_mg",model.c.clint_mic_ul_min_mg*2),
               ("ka_half","ka_h",model.c.ka_h*.5),("ka_double","ka_h",model.c.ka_h*2),
               ("Fa_half","fa",model.c.fa*.5),("all_Kp_half","kp",.5),("all_Kp_double","kp",2.)]
    for label,key,value in scenarios:
        c=copy.deepcopy(model.c)
        if key=="kp": c.kp_overrides={n:x*value for n,x in model.kp.items()}
        elif key: setattr(c,key,value)
        _,m,_=steady_state(PBPK(c,model.method),100,12)
        rows.append({"scenario":label,"status":m["status"],"coverage_pct":m.get("coverage_24h_pct"),
                     "cmax_mg_l":m.get("cmax_mg_l"),"ctrough_mg_l":m.get("trough_predose_mg_l"),
                     "dose_mg":100,"tau_h":12})
    return rows


def update_readme(out):
    path=out/"README.md"
    start,end="<!-- TASK4:BEGIN -->","<!-- TASK4:END -->"
    block=f"""{start}
# Developability Suite — Task 4

本目录是任务四的独立模块；默认结果均使用明确标注的示例输入。通过 JSON 模板接入指定先导化合物的参数，不自动把上游基准面板当作实测 PBPK 数据。

| Module | Input / output |
|---|---|
| 4A: IVIVE | Microsomal clearance and binding to whole-liver clearance. |
| 4B: Distribution | Seven amount states with a GI depot and serial pulmonary circulation. |
| 4C: Exposure | IV/oral profiles, AUC to infinity, multiple doses and periodic steady state. |
| 4D: Dose screening | QD/BID free-plasma coverage subject to a total-plasma peak constraint. |

## Run

Requires Python 3.10+, NumPy, SciPy and Matplotlib.

Run the commands below from the repository root. The existing calculation is archived in this project directory; use a separate output directory for a new calculation.

```sh
python -m pip install -r projects/task04_pbpk/requirements.txt
python projects/task04_pbpk/run_task4_pbpk_pharmacokinetics_dose_prediction.py --self-test --out work/task4_reproduction
python projects/task04_pbpk/run_task4_pbpk_pharmacokinetics_dose_prediction.py --write-example work/task4_compound.json
python projects/task04_pbpk/run_task4_pbpk_pharmacokinetics_dose_prediction.py --config work/task4_compound.json --out work/task4_custom
```

The JSON schema is the exported template: unknown keys and invalid values are rejected. Partial configurations override explicit defaults; omitted fields remain assumptions and are recorded in the resolved configuration. Set `km_unbound_mg_l` to a positive value for optional saturable hepatic metabolism. Nonlinear runs can take substantially longer and may have no finite steady state. Dose bounds are **per administration**; the default lattice is 1 mg. `target_tissue` defaults to `rest`, an explicit proxy. The toxicity threshold is **total plasma mg/L**; the efficacy threshold is **free plasma nM**.

## Deliverables

- [Self-contained script](run_task4_pbpk_pharmacokinetics_dose_prediction.py)
- [中文报告](PBPK_DOSE_PREDICTION_REPORT_ZH.md) / [English report](PBPK_DOSE_PREDICTION_REPORT_EN.md)
- [Resolved inputs](parameters_used.json) / [Editable example](compound_example.json)
- [PK and dose results](results_summary.json) / [Numerical checks](verification.json)
- `single_iv.csv`, `single_oral.csv`, `multidose_bid_7days.csv`, `steady_state_bid.csv`
- `dose_titration.csv`, `sensitivity_analysis.csv`, `manifest.json`
- `figures_task4/`: four requested publication-format PNGs at 300 DPI

The original brief mixes blood flow with plasma concentration and treats lung as a parallel organ. This implementation uses Rb consistently and places lung in series. The seven states include a GI **depot**, not a separate gut-wall tissue. Kp is a documented Poulin–Theil composition approximation with a logD ionization adaptation, not full Rodgers–Rowland. AUC includes its infinite tail; day 7 and true steady state are distinct. Clinical recommendation remains unset because the required validation and toxicology evidence are absent.

## Git delivery

To stage, commit and push generated outputs in an existing repository on local `main` with `origin` configured, run the script with `--git-sync`. It stages only this output directory and the script, refuses pre-existing staged changes, uses the requested commit message, and pushes `origin main` without force. Failures are saved in `git_delivery_status.json` and return a nonzero exit status. The ordinary run has no Git side effects.
{end}
"""
    old=path.read_text(encoding="utf-8") if path.exists() else ""
    if start in old and end in old:
        old=old[:old.index(start)]+block+old[old.index(end)+len(end):]
    else:
        old=old.rstrip()+"\n\n"+block if old else block
    path.write_text(old.rstrip()+"\n",encoding="utf-8", newline="\n")


def git_sync(out):
    def git(*args):
        return subprocess.check_output(["git","-C",str(out),*args],text=True,stderr=subprocess.STDOUT).strip()
    root=Path(git("rev-parse","--show-toplevel")).resolve()
    if git("branch","--show-current")!="main":
        raise RuntimeError("Git sync requires the existing local main branch")
    git("remote","get-url","origin")
    if git("diff","--cached","--name-only"):
        raise RuntimeError("Git sync found pre-existing staged changes; refusing to include unrelated work")
    paths=[str(out.relative_to(root))]
    script=Path(__file__).resolve()
    if script.is_relative_to(root) and not script.is_relative_to(out):
        paths.append(str(script.relative_to(root)))
    git("add","--",*paths)
    committed=False
    if git("diff","--cached","--name-only"):
        git("commit","-m",COMMIT_MESSAGE); committed=True
    git("push","origin","main")
    return {"status":"pushed","commit_created":committed,"commit":git("rev-parse","HEAD")}


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config",type=Path)
    parser.add_argument("--out",type=Path,default=Path.cwd())
    parser.add_argument("--write-example",type=Path)
    parser.add_argument("--method",choices=("Radau","BDF"),default="Radau")
    parser.add_argument("--self-test",action="store_true",help="Run independent regression checks before generating outputs")
    parser.add_argument("--git-sync",action="store_true",help="Stage generated output, commit and push an existing main/origin repository")
    args=parser.parse_args(argv)
    if args.write_example:
        args.write_example.parent.mkdir(parents=True,exist_ok=True)
        write_json(args.write_example,asdict(Config()))
        print(f"Wrote editable example: {args.write_example}")
        return 0
    data={}
    if args.config:
        data=json.loads(args.config.read_text(encoding="utf-8-sig"))
        if not isinstance(data,dict): raise ValueError("Config must be a JSON object")
        unknown=set(data)-set(Config.__dataclass_fields__)
        if unknown: raise ValueError(f"Unknown config keys: {sorted(unknown)}")
    c=Config(**data); c.validate()
    out=args.out.resolve(); out.mkdir(parents=True,exist_ok=True)
    print("Research PBPK run; all default compound and toxicity inputs are synthetic.",flush=True)
    checks=self_tests() if args.self_test else []
    model=PBPK(c,args.method)
    print("Integrating IV/oral doses and seven-day dosing...",flush=True)
    iv,im=single_dose(model,"iv"); oral,om=single_dose(model,"oral")
    mt,my,intervals=multidose(model)
    ss=steady_state(model,100,12)
    day,dayturns=extrema(intervals[-1],model,12)
    day["coverage_last_interval_pct"]=coverage(intervals[-1],model,12,turns=dayturns)
    day["coverage_last_24h_pct"]=(coverage(intervals[-2],model,12)+day["coverage_last_interval_pct"])/2
    day["within_1pct_of_ss"]=False
    day["state_error_vs_ss_relative"]=None
    if ss[1]["status"]=="converged":
        first,_=extrema(intervals[0],model,12)
        ss[1]["accumulation_ratio"]=(ss[1]["cmax_mg_l"]/first["cmax_mg_l"] if first["cmax_mg_l"] else None)
        reference=ss[0].y[:7,-1]
        state_error=float(np.max(abs(intervals[-1].y[:7,-1]-reference))/max(1e-12,np.max(abs(reference))))
        day["state_error_vs_ss_relative"]=state_error
        errors=[abs(day[k]-ss[1][k])/max(1e-12,ss[1][k]) for k in ("cmax_mg_l","trough_predose_mg_l")]
        day["within_1pct_of_ss"]=bool(max([state_error,*errors])<=.01)
    print("Screening QD/BID doses and parameter sensitivity...",flush=True)
    optimization,rows,curves=optimize_doses(model)
    sensitivities=sensitivity(model)
    result={"compound":c.compound_name,"synthetic_inputs":c.synthetic_inputs,
            "model_type":"linear" if model.linear else "saturable_hepatic_metabolism",
            "input_keys_provided":sorted(data),"input_keys_defaulted":sorted(set(asdict(c))-set(data)),
            "ivive":model.ivive(),"partition_info":model.partition_info,"kp_tissue_plasma":model.kp,
            "single_dose":{"iv":im,"oral":om},
            "oral_bioavailability_auc_ratio_pct":100*om["auc_0_inf_mg_h_l"]/im["auc_0_inf_mg_h_l"],
            "bioavailability_interpretation":"absolute modeled systemic fraction" if model.linear else "apparent equal-dose AUC ratio, not absolute bioavailability",
            "day7_last_interval":day,"bid_100mg_steady_state":ss[1],"dose_optimization":optimization}
    write_json(out/"parameters_used.json",asdict(c)); write_json(out/"compound_example.json",asdict(Config()))
    write_json(out/"results_summary.json",result)
    for name,sol in (("single_iv",iv),("single_oral",oral)):
        t=np.unique(np.r_[0,np.geomspace(1e-5,48,500),np.linspace(0,48,1441)])
        write_csv(out/f"{name}.csv",profile_rows(t,sol.sol(t),model))
    write_csv(out/"multidose_bid_7days.csv",profile_rows(mt,my,model))
    if ss[0] is not None:
        t=np.linspace(0,12,1441)
        write_csv(out/"steady_state_bid.csv",profile_rows(t,ss[0].sol(t),model))
    elif (out/"steady_state_bid.csv").exists():
        (out/"steady_state_bid.csv").unlink()  # prevent stale successful SS output
    write_csv(out/"dose_titration.csv",rows); write_csv(out/"sensitivity_analysis.csv",sensitivities)
    verification={"regression_tests_run":bool(args.self_test),"checks":checks,
                  "runtime_checks":["solver success","finite nonnegative amounts","mass conservation", "finite AUC tail", "periodic boundary when SS exists","dose feasibility and predecessor"],
                  "single_iv_max_mass_error_mg":im["mass_balance_max_error_mg"],
                  "single_oral_max_mass_error_mg":om["mass_balance_max_error_mg"],
                  "multidose_final_mass_error_mg":float(abs(my[:10,-1].sum()-1400)),
                  "clinical_validation_performed":False}
    write_json(out/"verification.json",verification)
    print("Writing bilingual reports and four 300-DPI figures...",flush=True)
    reports(out,model,result); update_readme(out)
    (out/"requirements.txt").write_text("numpy>=1.24\nscipy>=1.10\nmatplotlib>=3.7\n",encoding="utf-8", newline="\n")
    make_figures(out,model,iv,oral,result["single_dose"],mt,my,ss,curves)
    manifest={"python":platform.python_version(),"numpy":np.__version__,"scipy":scipy.__version__,
              "matplotlib":matplotlib.__version__,"solver":args.method,
              "rtol":c.rtol,"atol":c.atol,"script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "files_sha256":{}}
    for p in sorted(out.rglob("*")):
        if p.is_file() and ".git" not in p.parts and "__pycache__" not in p.parts and p.name not in {"manifest.json","git_delivery_status.json"}:
            manifest["files_sha256"][p.relative_to(out).as_posix()]=hashlib.sha256(p.read_bytes()).hexdigest()
    write_json(out/"manifest.json",manifest)
    if args.git_sync:
        try:
            status=git_sync(out)
        except (subprocess.CalledProcessError,FileNotFoundError,RuntimeError,ValueError) as exc:
            detail=getattr(exc,"output",None)
            status={"status":"failed","reason":str(exc),"detail":detail}
            write_json(out/"git_delivery_status.json",status)
            print(f"Artifacts generated; Git delivery failed: {exc}",file=sys.stderr)
            return 2
        write_json(out/"git_delivery_status.json",status)
    else:
        write_json(out/"git_delivery_status.json",{"status":"not_requested_for_this_run","note":"This run did not request Git changes. Use --git-sync in an existing repository when publication is intended."})
    print(json.dumps({"output_directory":str(out),"tests_passed":len(checks),
                      "selected_research_regimen":optimization["selected_research_regimen"]},ensure_ascii=False,indent=2),flush=True)
    return 0


if __name__=="__main__":
    try:
        raise SystemExit(main())
    except (ValueError,TypeError,RuntimeError,AssertionError,np.linalg.LinAlgError) as exc:
        print(f"ERROR: {exc}",file=sys.stderr)
        raise SystemExit(1)
