#!/usr/bin/env python3
"""Task 10: reproducible RNA ensembles and hypothetical linked equilibria.

NMR coordinates are experimental; derived geometry is a descriptor. Binding
energies and the mapping of U1 occupancy to inclusion are assumptions, not a
prediction of clinical risdiplam dose, splice efficacy, or a kinetic mechanism.
"""
from __future__ import annotations
import argparse
import copy
import csv
import hashlib
import json
import math
import platform
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
import scipy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
R = 0.00198720425864083  # kcal mol-1 K-1
COLORS = ["#167d9a", "#df8837", "#89949d", "#8753a4"]
FIGURES = ["fig1_rna_secondary_basepair_entropy_map.png",
           "fig2_rna_3d_groove_electrostatic_cleft.png",
           "fig3_coupled_thermodynamic_binding_cycle.png",
           "fig4_alternative_splicing_exon_inclusion.png"]


def defaults():
    return {"temperature_C": 37.0, "salt_M": 0.15,
            "sequence": "AUACUUACCUGUUCGGGAGUAAGUCU", "mutant_position_1based": 18,
            "mutant_base": "C", "temperature_scan_C": [20, 25, 30, 37, 42, 50, 60],
            "sequence_status": "Engineered tethered hairpin: 6HMI A + UUCG + B; PSU replaced by U. Not native SMN2 pre-mRNA.",
            "rna_target_total_M": 2e-9, "rna_offtarget_total_M": 3e-8,
            "u1_total_M": 1e-7, "u1_kd_M": 2e-6,
            "dose_min_M": 1e-10, "dose_max_M": 1e-2, "dose_points": 121,
            "ligands": [
                {"name": "Risdiplam-like model", "conf_kcal": 5.5, "stack": -7., "hbond": -4., "electrostatic": -3., "desolvation": 1., "entropy": 1., "closed_kd_M": .05, "splice_kcal": -3.5, "offtarget_kd_factor": 100., "offtarget_splice_kcal": -1.},
                {"name": "Aminoglycoside-like model", "conf_kcal": 4., "stack": -3., "hbond": -3., "electrostatic": -8., "desolvation": 3., "entropy": 2., "closed_kd_M": .005, "splice_kcal": -2., "offtarget_kd_factor": 2., "offtarget_splice_kcal": -2.},
                {"name": "Inactive mismatch model", "conf_kcal": 5.5, "stack": -2., "hbond": -1., "electrostatic": -1., "desolvation": 1., "entropy": 1., "closed_kd_M": 1., "splice_kcal": 0., "offtarget_kd_factor": 1., "offtarget_splice_kcal": 0.},
                {"name": "Flexible bis-intercalator model", "conf_kcal": 6.5, "stack": -9., "hbond": -4., "electrostatic": -2., "desolvation": 1.5, "entropy": 3.5, "closed_kd_M": .01, "splice_kcal": -3., "offtarget_kd_factor": 5., "offtarget_splice_kcal": -2.5}
            ]}


def validate_config(c):
    if set(c) != set(defaults()):
        raise ValueError("Configuration keys must match --write-example exactly")
    s = c["sequence"]
    if not isinstance(s, str) or not 8 <= len(s) <= 150 or set(s) - set("ACGU"):
        raise ValueError("sequence must contain 8–150 canonical A/C/G/U bases")
    if not isinstance(c["mutant_position_1based"], int) or not 1 <= c["mutant_position_1based"] <= len(s) or c["mutant_base"] not in {"A", "C", "G", "U"}:
        raise ValueError("Invalid mutant position or base")
    for k in ("temperature_C",):
        if not math.isfinite(c[k]) or not 0 <= c[k] <= 90:
            raise ValueError("Temperature outside 0–90 C")
    if not c["temperature_scan_C"] or any(not 0 <= t <= 90 for t in c["temperature_scan_C"]):
        raise ValueError("Invalid temperature scan")
    for k in ("salt_M", "rna_target_total_M", "u1_kd_M", "dose_min_M", "dose_max_M"):
        if not math.isfinite(c[k]) or c[k] <= 0:
            raise ValueError(k + " must be positive and finite")
    for k in ("rna_offtarget_total_M", "u1_total_M"):
        if not math.isfinite(c[k]) or c[k] < 0:
            raise ValueError(k + " must be nonnegative and finite")
    if c["dose_max_M"] <= c["dose_min_M"] or not isinstance(c["dose_points"], int) or not 20 <= c["dose_points"] <= 1001:
        raise ValueError("Invalid dose grid")
    if len(c["ligands"]) != 4 or len({x["name"] for x in c["ligands"]}) != 4:
        raise ValueError("Exactly four uniquely named ligand models required")
    keys = set(defaults()["ligands"][0])
    for x in c["ligands"]:
        if set(x) != keys or any(not math.isfinite(x[k]) for k in keys - {"name"}):
            raise ValueError("Invalid ligand fields")
        if x["closed_kd_M"] <= 0 or x["offtarget_kd_factor"] <= 0 or not 0 <= x["conf_kcal"] <= 20:
            raise ValueError("Invalid ligand affinity or conformational penalty")
        if any(abs(x[k]) > 30 for k in ["stack", "hbond", "electrostatic", "desolvation", "entropy", "splice_kcal", "offtarget_splice_kcal"]):
            raise ValueError("Free energy outside bounded model domain")
    return c


def save_json(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False, default=lambda x: x.item() if isinstance(x, np.generic) else str(x)) + "\n", encoding="utf-8")


def save_csv(path, rows):
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def rna_module():
    try:
        import RNA
    except ImportError as e:
        raise RuntimeError("ViennaRNA required; install requirements or pass --rna-library PATH. No toy fallback is used.") from e
    return RNA


def secondary_ensemble(sequence, temperature_C, salt_M=.15, dangles=2):
    RNA = rna_module()
    md = RNA.md()
    md.temperature = temperature_C
    md.salt = salt_M
    md.dangles = dangles
    md.compute_bpp = 1
    fc = RNA.fold_compound(sequence, md)
    mfe_structure, mfe = fc.mfe()
    fc.exp_params_rescale(mfe)
    _, free_energy = fc.pf()
    upper = np.asarray(fc.bpp(), dtype=float)[1:, 1:]
    p = upper + upper.T
    unpaired = 1. - p.sum(axis=1)
    if unpaired.min() < -1e-8:
        raise ArithmeticError("RNA pairing probabilities fail normalization")
    unpaired = np.maximum(unpaired, 0.)
    states = np.column_stack((p, unpaired))
    ent = -(states * np.log2(np.maximum(states, 1e-300))).sum(axis=1)
    return {"sequence": sequence, "temperature_C": temperature_C,
            "mfe_structure": mfe_structure, "mfe_kcal_mol": float(mfe),
            "ensemble_free_energy_kcal_mol": float(free_energy),
            "partition_kT_kcal_mol": float(fc.exp_params.kT/1000.),
            "log_partition": float(-free_energy / (fc.exp_params.kT/1000.)),
            "p": p, "unpaired": unpaired, "entropy": ent}


def enumerate_structures(sequence):
    """Exhaustive noncrossing canonical/wobble structures, hairpin turn >=3."""
    pairs = {"AU", "UA", "CG", "GC", "GU", "UG"}
    def rec(lo, hi):
        if lo > hi:
            return [""]
        ans = ["." + v for v in rec(lo + 1, hi)]
        for j in range(lo + 4, hi + 1):
            if sequence[lo] + sequence[j] in pairs:
                ans.extend("(" + a + ")" + b for a in rec(lo + 1, j - 1) for b in rec(j + 1, hi))
        return ans
    return rec(0, len(sequence)-1)


def read_pdb(path):
    models, current = [], []
    for line in path.read_text(encoding="ascii").splitlines():
        if line.startswith("ENDMDL"):
            models.append(current)
            current = []
        elif line.startswith(("ATOM  ", "HETATM")) and line[16:17] in (" ", "A"):
            element = line[76:78].strip() or line[12:16].strip()[0]
            if element != "H":
                current.append({"atom": line[12:16].strip(), "resname": line[17:20].strip(),
                                "chain": line[21], "resid": int(line[22:26]), "element": element,
                                "xyz": np.array([float(line[a:a+8]) for a in (30, 38, 46)])})
    if current:
        models.append(current)
    if not models:
        raise ValueError("No PDB coordinates")
    return models


def debye_length_angstrom(salt_M, temperature_K, dielectric=78.5):
    # SI expression: lambda_D = sqrt(eps eps0 kBT/(2 NA 1000 I e^2)).
    return math.sqrt(dielectric * 8.8541878128e-12 * 1.380649e-23 * temperature_K /
                     (2 * 6.02214076e23 * 1000 * salt_M * 1.602176634e-19 ** 2)) * 1e10


def phosphate_potential_mV(point, phosphates, salt_M, temperature_K):
    d = np.linalg.norm(phosphates-point, axis=1)
    if np.any(d < .01):
        raise ValueError("Electrostatic probe overlaps a point charge")
    lam = debye_length_angstrom(salt_M, temperature_K)
    return float(-14399.6454784255 / 78.5 * np.sum(np.exp(-d/lam)/d))


EDGES = {
    "A": {"Watson-Crick": {"N1": "acceptor", "N6": "donor"}, "Hoogsteen": {"N7": "acceptor", "N6": "donor"}, "Sugar": {"N3": "acceptor"}},
    "G": {"Watson-Crick": {"O6": "acceptor", "N1": "donor", "N2": "donor"}, "Hoogsteen": {"N7": "acceptor", "O6": "acceptor"}, "Sugar": {"N3": "acceptor", "N2": "donor"}},
    "C": {"Watson-Crick": {"N3": "acceptor", "N4": "donor", "O2": "acceptor"}, "Hoogsteen": {"N4": "donor"}, "Sugar": {"O2": "acceptor"}},
    "U": {"Watson-Crick": {"O4": "acceptor", "N3": "donor", "O2": "acceptor"}, "Hoogsteen": {"O4": "acceptor"}, "Sugar": {"O2": "acceptor"}},
    "PSU": {"Watson-Crick": {"O4": "acceptor", "N3": "donor", "O2": "acceptor"}, "Hoogsteen": {"O4": "acceptor"}, "Sugar": {"O2": "acceptor", "N1": "donor"}}}
VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "P": 1.80, "F": 1.47}
PAIR_MAP = [(3, 20), (4, 19), (5, 18), (6, 17), (7, 16), (8, 15), (9, 13)]  # PDB author residue IDs


def geometry(c, source_dir):
    """Observed NMR ensembles; fixed operational cleft descriptors, not Curves+."""
    rows, annotation, stackrows, render = [], [], [], {}
    for pdbid in ("6HMI", "6HMO"):
        models = read_pdb(source_dir / (pdbid + ".pdb"))
        render[pdbid] = models[0]
        for im, atoms in enumerate(models, 1):
            rna = [a for a in atoms if a["resname"] in EDGES]
            coords = np.array([a["xyz"] for a in rna])
            radii = np.array([VDW.get(a["element"], 1.7) for a in rna])
            phosphates = np.array([a["xyz"] for a in rna if a["atom"] == "P"])
            residues = {}
            for a in rna:
                residues.setdefault((a["chain"], a["resid"]), []).append(a)
                if im == 1:
                    for edge, sites in EDGES[a["resname"]].items():
                        if a["atom"] in sites:
                            annotation.append({"pdb": pdbid, "chain": a["chain"], "resid": a["resid"], "base": a["resname"], "atom": a["atom"], "face": edge, "role": sites[a["atom"]], "x_A": a["xyz"][0], "y_A": a["xyz"][1], "z_A": a["xyz"][2]})
            centers = []
            for ai, bi in PAIR_MAP:
                both = residues[("A", ai)] + residues[("B", bi)]
                base = [a["xyz"] for a in both if "'" not in a["atom"] and a["atom"] not in ("P", "OP1", "OP2")]
                centers.append(np.mean(base, axis=0))
            centers = np.array(centers)
            _, _, vt = np.linalg.svd(centers-centers.mean(axis=0), full_matrices=False)
            axis = vt[0]
            if np.dot(axis, centers[-1]-centers[0]) < 0:
                axis = -axis
            for station, (ai, bi) in enumerate(PAIR_MAP):
                both = residues[("A", ai)] + residues[("B", bi)]
                center = centers[station]
                for face, label in [("Hoogsteen", "major-facing"), ("Sugar", "minor-facing")]:
                    face_xyz = [a["xyz"] for a in both if a["atom"] in EDGES[a["resname"]][face]]
                    edgecenter = np.mean(face_xyz, axis=0)
                    direction = edgecenter - center
                    direction -= np.dot(direction, axis) * axis
                    direction /= np.linalg.norm(direction)
                    probe = edgecenter + 4. * direction
                    closest = []
                    for chain in ("A", "B"):
                        pp = np.array([a["xyz"] for a in rna if a["atom"] == "P" and a["chain"] == chain])
                        closest.append(pp[np.argmin(np.linalg.norm(pp-probe, axis=1))])
                    ppdist = np.linalg.norm(closest[0]-closest[1])
                    clearance = float(np.min(np.linalg.norm(coords-probe, axis=1)-radii))
                    # Diameter of largest sphere at a specified probe point. This is
                    # distinct from the P–P chord and from a standard groove width.
                    rows.append({"pdb": pdbid, "model": im, "station": station+1,
                                 "A_resid": ai, "B_resid": bi, "face": label,
                                 "axis_A": float(np.dot(center-centers[0], axis)),
                                 "phosphate_chord_A": float(ppdist),
                                 "vdw_sphere_diameter_A": 2*max(0., clearance),
                                 "water_probe_clearance_A": clearance-1.4,
                                 "edge_to_backbone_depth_A": float(abs(np.dot(np.mean(closest, axis=0)-edgecenter, direction))),
                                 "DH_potential_mV": phosphate_potential_mV(probe, phosphates, c["salt_M"], c["temperature_C"]+273.15)})
            # Bulged B14 adenine (third base of chain B) and purine neighbours.
            rings = {}
            for key, aa in residues.items():
                base = aa[0]["resname"]
                if base in ("A", "G"):
                    rr = np.array([a["xyz"] for a in aa if a["atom"] in {"N9", "C8", "N7", "C5", "C6", "N1", "C2", "N3", "C4"}])
                    centroid = rr.mean(axis=0)
                    _, _, v = np.linalg.svd(rr-centroid, full_matrices=False)
                    rings[key] = (centroid, v[-1])
            for other in [("B", 13), ("B", 15), ("A", 1), ("A", 3), ("A", 7), ("A", 11)]:
                if other not in rings or ("B", 14) not in rings:
                    continue
                a, n = rings[("B", 14)]
                b, m = rings[other]
                stackrows.append({"pdb": pdbid, "model": im, "base1": "B14:A", "base2": f"{other[0]}{other[1]}", "centroid_distance_A": float(np.linalg.norm(a-b)), "normal_angle_deg": float(np.degrees(np.arccos(np.clip(abs(np.dot(n, m)), 0, 1))))})
    return rows, annotation, stackrows, render


def ligand_parameters(x, temperature_K):
    rt = R*temperature_K
    intrinsic = sum(x[k] for k in ("stack", "hbond", "electrostatic", "desolvation", "entropy"))
    q = math.exp(-x["conf_kcal"]/rt)
    ko = math.exp(intrinsic/rt)  # M, using 1 M standard state
    kc = x["closed_kd_M"]
    pc = q/(1+q)
    keff = (1+q)/(1/kc+q/ko)
    return {"name": x["name"], "conf_kcal": x["conf_kcal"], "intrinsic_bind_kcal": intrinsic,
            "competent_fraction": pc, "intrinsic_kd_M": ko,
            "open_only_kdeff_M": ko/pc, "full_kdeff_M": keff,
            "net_standard_bind_kcal": rt*math.log(keff),
            "closed_bind_kcal": rt*math.log(kc),
            "bound_conf_kcal": x["conf_kcal"]+intrinsic-rt*math.log(kc),
            "cycle_residual_kcal": (x["conf_kcal"]+intrinsic)-(rt*math.log(kc)+(x["conf_kcal"]+intrinsic-rt*math.log(kc)))}


def state_weights(free_l, free_u, q, ko, kc, ku, alpha):
    # [C, O, CL, OL, CU, OU, CLU, OLU]. Conditional reciprocal
    # affinity shifts share alpha, which enforces microscopic reversibility.
    l_c, l_o, u = free_l/kc, free_l/ko, free_u/ku
    w = np.array([1., q, l_c, q*l_o, u, q*u, l_c*u, q*l_o*u*alpha])
    return w/w.sum()


def make_sites(c, ligand, mutant=False):
    rt = R*(c["temperature_C"]+273.15)
    p = ligand_parameters(ligand, c["temperature_C"]+273.15)
    q = math.exp(-(ligand["conf_kcal"] + (1.5 if mutant else 0.))/rt)
    target = (c["rna_target_total_M"], q, p["intrinsic_kd_M"]*(10 if mutant else 1),
              ligand["closed_kd_M"], c["u1_kd_M"]*(3 if mutant else 1),
              math.exp(-(max(ligand["splice_kcal"], -1.) if mutant else ligand["splice_kcal"])/rt))
    off = (c["rna_offtarget_total_M"], math.exp(-(max(0., ligand["conf_kcal"]-1.5))/rt),
           p["intrinsic_kd_M"]*ligand["offtarget_kd_factor"], ligand["closed_kd_M"],
           c["u1_kd_M"]*2., math.exp(-ligand["offtarget_splice_kcal"]/rt))
    return [target, off]


def solve_equilibrium(total_l, total_u, sites):
    def probs(l, u):
        return [state_weights(l, u, *s[1:]) for s in sites]
    def bound(p, indices):
        return sum(s[0]*z[indices].sum() for s, z in zip(sites, p))
    def solve_u(l):
        if total_u == 0:
            return 0.
        return brentq(lambda u: u+bound(probs(l, u), [4, 5, 6, 7])-total_u,
                      0., total_u, xtol=1e-23, rtol=1e-12)
    if total_l == 0:
        l = 0.
    else:
        l = brentq(lambda ll: ll+bound(probs(ll, solve_u(ll)), [2, 3, 6, 7])-total_l,
                   0., total_l, xtol=1e-23, rtol=1e-12)
    u = solve_u(l)
    p = probs(l, u)
    return {"free_l_M": l, "free_u_M": u, "p": p,
            "ligand_mass_error_M": l+bound(p, [2, 3, 6, 7])-total_l,
            "u1_mass_error_M": u+bound(p, [4, 5, 6, 7])-total_u}


def saturation_probabilities(total_u, sites):
    """Exact L -> infinity state weights, with finite U1 mass conservation."""
    def probs(u):
        ans = []
        for _, q, ko, kc, ku, alpha in sites:
            a, b = 1/kc, q/ko
            w = np.array([0., 0., a, b, 0., 0., a*u/ku, b*u*alpha/ku])
            ans.append(w/w.sum())
        return ans
    u = 0. if total_u == 0 else brentq(lambda uf: uf+sum(s[0]*p[4:].sum() for s, p in zip(sites, probs(uf)))-total_u,
                                      0., total_u, xtol=1e-23, rtol=1e-12)
    return probs(u)


def contiguous_grid_bands(rows, baseline_off):
    bands, current = [], []
    for row in rows:
        eligible = row["inclusion_proxy_percent"] >= 85 and row["offtarget_inclusion_proxy_percent"]-baseline_off <= 10
        if eligible:
            current.append(row["total_l_M"])
        elif current:
            bands.append({"min_M": current[0], "max_M": current[-1], "grid_points": len(current), "right_censored_at_grid_max": False})
            current = []
    if current:
        bands.append({"min_M": current[0], "max_M": current[-1], "grid_points": len(current), "right_censored_at_grid_max": True})
    return bands


def splicing_scan(c):
    doses = np.r_[0., np.geomspace(c["dose_min_M"], c["dose_max_M"], c["dose_points"])]
    rows, endpoints = [], []
    for x in c["ligands"]:
        for mutant in (False, True):
            sites = make_sites(c, x, mutant)
            baseline = solve_equilibrium(0, c["u1_total_M"], sites)["p"][0][4:].sum()
            local = []
            for dose in doses:
                ans = solve_equilibrium(float(dose), c["u1_total_M"], sites)
                pt, po = ans["p"]
                row = {"ligand": x["name"], "site": "mutant scenario" if mutant else "WT scenario", "total_l_M": float(dose),
                       "free_l_M": ans["free_l_M"], "free_u_M": ans["free_u_M"],
                       "target_U1_occupancy": float(pt[4:].sum()), "inclusion_proxy_percent": float(100*pt[4:].sum()),
                       "offtarget_inclusion_proxy_percent": float(100*po[4:].sum()),
                       "target_ligand_occupancy": float(pt[[2, 3, 6, 7]].sum()),
                       "ligand_mass_error_M": ans["ligand_mass_error_M"], "u1_mass_error_M": ans["u1_mass_error_M"]}
                rows.append(row)
                local.append(row)
            asym = saturation_probabilities(c["u1_total_M"], sites)[0][4:].sum()
            midpoint = (baseline+asym)/2
            ec50, crossings = None, []
            differences = np.diff([a["target_U1_occupancy"] for a in local])
            direction = ("flat" if np.max(np.abs(differences)) < 1e-9 else
                         "monotonic_increasing" if differences.min() >= -1e-9 else
                         "monotonic_decreasing" if differences.max() <= 1e-9 else "nonmonotonic")
            if abs(asym-baseline) > 1e-6:
                for lo, hi in zip(local[:-1], local[1:]):
                    if (lo["target_U1_occupancy"]-midpoint)*(hi["target_U1_occupancy"]-midpoint) <= 0:
                        crossing = brentq(lambda dl: solve_equilibrium(dl, c["u1_total_M"], sites)["p"][0][4:].sum()-midpoint,
                                          lo["total_l_M"], hi["total_l_M"], xtol=1e-19, rtol=1e-10)
                        crossings.append({"total_l_M": crossing, "direction": "up" if hi["target_U1_occupancy"] > lo["target_U1_occupancy"] else "down"})
            if direction == "monotonic_increasing" and asym > baseline+1e-6 and len(crossings) == 1:
                ec50 = crossings[0]["total_l_M"]
            baseline_off = local[0]["offtarget_inclusion_proxy_percent"]
            endpoints.append({"ligand": x["name"], "site": local[0]["site"], "baseline_proxy_percent": float(100*baseline),
                              "asymptotic_proxy_percent": float(100*asym), "grid_max_proxy_percent": max(a["inclusion_proxy_percent"] for a in local),
                              "ec50_total_M": ec50, "ec50_status": "estimated" if ec50 is not None else ("not_applicable_"+direction if direction != "monotonic_increasing" else "outside_grid"),
                              "response_direction_on_grid": direction, "half_dynamic_range_crossings": crossings,
                              "illustrative_selectivity_grid_bands": contiguous_grid_bands(local, baseline_off),
                              "band_status": "Contiguous qualifying sampled points only; no interpolation or clinical window"})
    return rows, endpoints


def numerical_checks(c):
    checks = {}
    RNA = rna_module()
    seq = "GCGAAAGC"
    md = RNA.md()
    md.temperature = c["temperature_C"]
    md.salt = c["salt_M"]
    md.dangles = 0
    fc = RNA.fold_compound(seq, md)
    structs = enumerate_structures(seq)
    energies = [fc.eval_structure(s) for s in structs]
    calc = secondary_ensemble(seq, c["temperature_C"], c["salt_M"], dangles=0)
    z = sum(math.exp(-e/calc["partition_kT_kcal_mol"]) for e in energies)
    checks["enumerated_structure_count"] = len(structs)
    checks["partition_log_absolute_error"] = abs(math.log(z)-calc["log_partition"])
    checks["partition_enumeration_pass"] = checks["partition_log_absolute_error"] < 1e-5
    bare = secondary_ensemble("AAAAAAAA", 37, c["salt_M"])
    checks["unpairable_sequence_pass"] = bool(np.max(bare["p"]) == 0 and np.max(bare["entropy"]) == 0)
    sites = make_sites(c, c["ligands"][0])
    e = solve_equilibrium(1e-8, c["u1_total_M"], sites)
    checks["finite_pool_mass_balance_pass"] = abs(e["ligand_mass_error_M"]) < 1e-16 and abs(e["u1_mass_error_M"]) < 1e-16
    checks["ligand_depletion_pass"] = e["free_l_M"] <= 1e-8
    no_u = solve_equilibrium(1e-5, 0, sites)
    checks["zero_u1_pass"] = bool(no_u["p"][0][4:].sum() == 0)
    p = state_weights(1e-6, 1e-7, .01, 1e-8, 1e-3, 1e-6, 100.)
    checks["thermodynamic_cycle_ratio_error"] = float(abs((p[7]*p[1])/(p[3]*p[5])-100.))
    checks["thermodynamic_cycle_pass"] = checks["thermodynamic_cycle_ratio_error"] < 1e-10
    plain = state_weights(0, 1e-7, .01, 1e-8, 1e-3, 1e-6, 100.)
    checks["no_ligand_baseline_pass"] = abs(plain[4:].sum()-1/11) < 1e-12
    checks["all_passed"] = all(v for k, v in checks.items() if k.endswith("_pass"))
    if not checks["all_passed"]:
        raise AssertionError(checks)
    return {k: v.item() if isinstance(v, np.generic) else v for k, v in checks.items()}


def figures(out, ens, geometry_rows, render, binding, doses, c):
    fd = out/"figures_task10"
    fd.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False, "savefig.dpi": 300})
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), gridspec_kw={"width_ratios": [1, 1.25]}, layout="constrained")
    im = axes[0].imshow(ens["p"], origin="lower", vmin=0, vmax=1, cmap="viridis", extent=(.5, len(ens["sequence"])+.5, .5, len(ens["sequence"])+.5))
    fig.colorbar(im, ax=axes[0], label="Pair probability")
    axes[0].set(xlabel="Nucleotide i", ylabel="Nucleotide j", title="ViennaRNA nearest-neighbour ensemble")
    pos = np.arange(1, len(ens["sequence"])+1)
    axes[1].bar(pos, ens["entropy"], color=COLORS[0], label="H: paired + unpaired states")
    axes[1].axhline(1.2, color=COLORS[1], ls="--", label="Requested 1.2-bit threshold")
    axes[1].set(xlabel="Engineered RNA position", ylabel="Shannon entropy (bits)", title="Unpaired is one state; entropy is not a flip barrier")
    axes[1].legend(fontsize=8)
    fig.suptitle(f"10A | Canonical-base tethered RNA model, {c['temperature_C']:g} C")
    fig.savefig(fd/FIGURES[0], dpi=300)
    plt.close(fig)
    fig = plt.figure(figsize=(13, 8), layout="constrained")
    ax = fig.add_subplot(221, projection="3d")
    for chain, color in zip(("A", "B"), COLORS):
        a = [a for a in render["6HMO"] if a["chain"] == chain and a["atom"] == "P"]
        xyz = np.array([v["xyz"] for v in a])
        ax.plot(*xyz.T, "o-", color=color, label=f"RNA {chain}", ms=3)
    lig = np.array([a["xyz"] for a in render["6HMO"] if a["resname"] == "GDZ"])
    ax.scatter(*lig.T, s=13, color=COLORS[3], label="SMN-C5 heavy atoms")
    ax.set(xlabel="x (A)", ylabel="y (A)", zlabel="z (A)", title="6HMO model 1: experimental coordinates")
    ax.legend(fontsize=8)
    for panel, key, ylabel in [(222, "phosphate_chord_A", "Nearest cross-strand P-P chord (A)"), (223, "vdw_sphere_diameter_A", "Local largest-sphere diameter (A)"), (224, "DH_potential_mV", "Screened phosphate potential (mV)")]:
        ax = fig.add_subplot(panel)
        for j, (pid, face) in enumerate([(p, f) for p in ("6HMI", "6HMO") for f in ("major-facing", "minor-facing")]):
            values = np.array([[r[key] for r in geometry_rows if r["pdb"] == pid and r["face"] == face and r["station"] == st] for st in range(1, 8)])
            mean, sd = values.mean(axis=1), values.std(axis=1)
            ax.plot(range(1, 8), mean, color=COLORS[j], label=pid+" "+face)
            lower = np.maximum(0., mean-sd) if key == "vdw_sphere_diameter_A" else mean-sd
            ax.fill_between(range(1, 8), lower, mean+sd, color=COLORS[j], alpha=.10)
        ax.set(xlabel="Paired station along duplex", ylabel=ylabel)
        if panel == 222:
            ax.legend(fontsize=7)
    fig.suptitle("10B | 20 NMR conformers per structure; spread is not a Boltzmann population\nGroove-facing descriptors are not standard groove widths or a PB solution")
    fig.savefig(fd/FIGURES[1], dpi=300)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), layout="constrained")
    xpos = np.arange(4)
    for offset, key, color, label in [(-.25, "conf_kcal", COLORS[1], "Opening penalty (assumed)"), (0, "intrinsic_bind_kcal", COLORS[0], "Intrinsic binding (assumed)"), (.25, "net_standard_bind_kcal", COLORS[3], "Full effective binding (calculated)")]:
        axes[0].bar(xpos+offset, [b[key] for b in binding], .23, color=color, label=label)
    axes[0].axhline(0, color="grey", lw=.6)
    axes[0].set(ylabel="Standard free energy (kcal/mol)", xticks=xpos, xticklabels=["Planar", "Polycation", "Mismatch", "Flexible"])
    axes[0].legend(fontsize=8)
    for i, b in enumerate(binding):
        axes[1].plot([i, i], [b["intrinsic_kd_M"]*1e6, b["full_kdeff_M"]*1e6], color=COLORS[i], lw=3)
        axes[1].scatter(i, b["intrinsic_kd_M"]*1e6, marker="o", color=COLORS[i])
        axes[1].scatter(i, b["full_kdeff_M"]*1e6, marker="s", color=COLORS[i])
    axes[1].set(yscale="log", ylabel="KD (uM): circle intrinsic; square effective", xticks=xpos, xticklabels=["Planar", "Polycation", "Mismatch", "Flexible"])
    fig.suptitle("10C | Four hypothetical scaffolds, 1 M standard state\nEquilibrium cycle does not distinguish induced fit from conformational selection")
    fig.savefig(fd/FIGURES[2], dpi=300)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), layout="constrained")
    for i, x in enumerate(c["ligands"]):
        for site, style in [("WT scenario", "-"), ("mutant scenario", "--")]:
            d = [a for a in doses if a["ligand"] == x["name"] and a["site"] == site and a["total_l_M"] > 0]
            axes[0].semilogx([a["total_l_M"]*1e6 for a in d], [a["inclusion_proxy_percent"] for a in d], color=COLORS[i], ls=style, label=x["name"].replace(" model", "") + (" WT" if style == "-" else " mutant"))
            if style == "-":
                axes[1].semilogx([a["total_l_M"]*1e6 for a in d], [a["offtarget_inclusion_proxy_percent"] for a in d], color=COLORS[i], label=x["name"].replace(" model", ""))
    for ax in axes:
        ax.set(xlabel="Total model ligand (uM), not administered dose", ylabel="U1 occupancy / inclusion proxy (%)", ylim=(0, 100))
    axes[0].axhline(85, ls=":", color="grey", label="85% target (not imposed)")
    axes[0].legend(fontsize=7, loc="upper left")
    axes[1].legend(fontsize=8)
    axes[0].set_title("Target: solid WT, dashed mutant")
    axes[1].set_title("One assumed competing off-target RNA pool")
    fig.suptitle("10D | Coupled finite-pool mass action; eight states per RNA\nOccupancy is not measured exon inclusion; no clinical window inferred")
    fig.savefig(fd/FIGURES[3], dpi=300)
    plt.close(fig)


def run(c, out, overwrite=False):
    started = time.perf_counter()
    if out.exists() and any(out.iterdir()) and not overwrite:
        raise FileExistsError("Output directory is nonempty; choose a new --out or explicitly pass --overwrite")
    out.mkdir(parents=True, exist_ok=True)
    save_json(out/"inputs"/"config.json", c)
    checks = numerical_checks(c)
    sequence = c["sequence"]
    idx = c["mutant_position_1based"]-1
    mutant = sequence[:idx]+c["mutant_base"]+sequence[idx+1:]
    secrows, secsummary = [], []
    ref = None
    for name, seq in [("WT_model", sequence), ("mutant_model", mutant)]:
        for temp in sorted(set(c["temperature_scan_C"]+[c["temperature_C"]])):
            e = secondary_ensemble(seq, temp, c["salt_M"])
            np.savetxt(out/(f"bpp_{name}_{temp:g}C.csv"), e["p"], delimiter=",", fmt="%.12g")
            secsummary.append({k: v for k, v in e.items() if k not in ("p", "unpaired", "entropy")} | {"variant": name, "flexible_positions_above_1p2_bits": [int(i+1) for i, v in enumerate(e["entropy"]) if v > 1.2]})
            secrows.extend({"variant": name, "temperature_C": temp, "position": i+1, "base": seq[i], "p_unpaired": e["unpaired"][i], "entropy_bits": e["entropy"][i]} for i in range(len(seq)))
            if name == "WT_model" and temp == c["temperature_C"]:
                ref = e
    save_csv(out/"secondary_position_statistics.csv", secrows)
    save_json(out/"secondary_ensemble_summary.json", secsummary)
    grows, edges, stacks, render = geometry(c, HERE/"inputs")
    save_csv(out/"groove_geometry.csv", grows)
    save_csv(out/"base_edge_annotations.csv", edges)
    save_csv(out/"bulge_stacking_geometry.csv", stacks)
    binding = [ligand_parameters(x, c["temperature_C"]+273.15) for x in c["ligands"]]
    save_csv(out/"thermodynamic_cycles.csv", binding)
    doses, endpoints = splicing_scan(c)
    save_csv(out/"splicing_mass_action.csv", doses)
    save_json(out/"splicing_endpoints.json", endpoints)
    sensitivity = []
    for conf in np.linspace(4, 7.5, 8):
        for delta in [-2., -2.5, -3., -3.5, -4.]:
            x = copy.deepcopy(c["ligands"][0])
            x["conf_kcal"], x["splice_kcal"] = float(conf), delta
            sites = make_sites(c, x)
            for dose in (1e-7, 1e-6, 1e-5):
                e = solve_equilibrium(dose, c["u1_total_M"], sites)
                sensitivity.append({"conf_kcal": conf, "splice_kcal": delta, "total_l_M": dose, "inclusion_proxy_percent": 100*float(e["p"][0][4:].sum())})
    save_csv(out/"assumption_sensitivity.csv", sensitivity)
    figures(out, ref, grows, render, binding, doses, c)
    checks["max_ligand_mass_error_M"] = max(abs(d["ligand_mass_error_M"]) for d in doses)
    checks["max_u1_mass_error_M"] = max(abs(d["u1_mass_error_M"]) for d in doses)
    checks["max_probability_normalization_error"] = float(np.max(abs(ref["p"].sum(axis=1)+ref["unpaired"]-1)))
    save_json(out/"verification.json", checks)
    summary = {"evidence_status": "Mixed public experimental coordinates, physical descriptors, and uncalibrated equilibrium scenarios. Not clinical predictions.",
               "rna_library": rna_module().__version__, "partition_runs": len(secsummary),
               "nmr_models_per_structure": {pid: len(read_pdb(HERE/"inputs"/(pid+".pdb"))) for pid in ("6HMI", "6HMO")},
               "geometry_rows": len(grows), "dose_response_rows": len(doses), "sensitivity_rows": len(sensitivity),
               "reference_mfe_structure": ref["mfe_structure"], "reference_ensemble_free_energy_kcal_mol": ref["ensemble_free_energy_kcal_mol"],
               "reference_max_entropy_bits": float(ref["entropy"].max()), "debye_length_A": debye_length_angstrom(c["salt_M"], c["temperature_C"]+273.15),
               "endpoints": endpoints, "verification": checks}
    save_json(out/"results_summary.json", summary)
    save_json(out/"run_log.json", {"seconds": time.perf_counter()-started, "python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__, "matplotlib": matplotlib.__version__, "ViennaRNA": rna_module().__version__, "source_script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
    freeze_manifest(out)
    return summary


def freeze_manifest(out):
    files = {p.relative_to(out).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(out.rglob("*")) if p.is_file() and p.name != "manifest.json" and "__pycache__" not in p.parts}
    save_json(out/"manifest.json", {"script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(), "files_sha256": files})


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=HERE/"results")
    ap.add_argument("--config", type=Path)
    ap.add_argument("--write-example", type=Path)
    ap.add_argument("--rna-library", type=Path, help="Optional isolated ViennaRNA package directory")
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--overwrite", action="store_true", help="Explicitly permit replacing generated outputs in --out")
    args = ap.parse_args()
    if args.rna_library:
        sys.path.insert(0, str(args.rna_library.resolve()))
    c = defaults() if args.config is None else json.loads(args.config.read_text(encoding="utf-8"))
    validate_config(c)
    if args.write_example:
        if args.write_example.exists() and not args.overwrite:
            raise FileExistsError("Example config already exists; choose a new path")
        save_json(args.write_example, c)
        return
    if args.self_test:
        print(json.dumps(numerical_checks(c), indent=2))
        return
    summary = run(c, args.out.resolve(), overwrite=args.overwrite)
    print(json.dumps({k: v for k, v in summary.items() if k not in ("endpoints", "verification")}, indent=2))


if __name__ == "__main__":
    main()
