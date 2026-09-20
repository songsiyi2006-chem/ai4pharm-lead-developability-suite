"""Read-only checks of the requested Task 8-10 deliverables and provenance.

Checks files, image metadata, source syntax, manifests and repository links.
Scientific invariants are tested separately by tests/test_task8.py through
tests/test_task10.py. Passing this tool is not biological validation.
"""
from pathlib import Path
import argparse
import ast
import json
import subprocess

from PIL import Image

from validate_repository_layout import (
    ROOT, PROJECTS, MANIFESTS, check_markdown_links, check_project_manifests,
    error, inventory, relative, sha, write_report,
)

SPECS = (
    (8, "ADC_TRANSLATIONAL_ENGINEERING_REPORT", (
        "fig1_adc_conjugation_dar_distribution.png",
        "fig2_analytical_hic_chromatogram_twin.png",
        "fig3_intracellular_lysosomal_release_ode.png",
        "fig4_bystander_killing_spatiotemporal_contour.png",
    )),
    (9, "CRYOEM_CRYPTIC_POCKET_REPORT", (
        "fig1_cryoem_latent_conformational_manifold.png",
        "fig2_msm_free_energy_pathway.png",
        "fig3_dynamic_pocket_volume_druggability.png",
        "fig4_prs_allosteric_network_matrix.png",
    )),
    (10, "RNA_TARGETED_CADD_REPORT", (
        "fig1_rna_secondary_basepair_entropy_map.png",
        "fig2_rna_3d_groove_electrostatic_cleft.png",
        "fig3_coupled_thermodynamic_binding_cycle.png",
        "fig4_alternative_splicing_exon_inclusion.png",
    )),
)


def inspect():
    errors = []
    result = {"validator": "validate_task8_10_delivery", "scope":
              "Requested files, source syntax, figure metadata, links and SHA256; not scientific validation"}
    try:
        _, included = inventory()
        records, figures = [], []
        for project, spec, manifest_spec in zip(PROJECTS[7:], SPECS, MANIFESTS[7:]):
            task, report, figure_names = spec
            folder = ROOT / "projects" / project
            driver = manifest_spec[3]
            required = ["README.md", "TASK.md", driver, "manifest.json",
                        f"{report}_ZH.md", f"{report}_EN.md"]
            required += [f"figures_task{task}/{name}" for name in figure_names]
            checked = []
            for name in required:
                path = folder / name
                if not path.is_file():
                    error(errors, "required_artifact_missing", path=relative(path))
                    continue
                if relative(path) not in included:
                    error(errors, "required_artifact_ignored", path=relative(path))
                checked.append({"path": relative(path), "sha256": sha(path)})
                if name.endswith(".py"):
                    ast.parse(path.read_text(encoding="utf-8-sig"), filename=relative(path))
            for name in figure_names:
                path = folder / f"figures_task{task}" / name
                if not path.is_file():
                    continue
                try:
                    with Image.open(path) as image:
                        dpi = image.info.get("dpi", ())
                        dimensions = list(image.size)
                        if len(dpi) != 2 or not all(abs(d - 300) < 0.02 for d in dpi):
                            error(errors, "figure_dpi", path=relative(path), actual=list(dpi))
                        image.verify()
                    figures.append({"path": relative(path), "sha256": sha(path),
                                    "dpi": list(dpi), "pixels": dimensions})
                except (OSError, ValueError) as exc:
                    error(errors, "invalid_figure", path=relative(path), detail=str(exc))
            data_files = sorted(name for name in included
                                if name.startswith(relative(folder) + "/")
                                and Path(name).suffix.lower() in {".csv", ".npz", ".npy", ".mrc"})
            if not data_files:
                error(errors, "numerical_data_missing", project=project)
            records.append({"project": project, "required_files": checked,
                            "numerical_data_files": data_files})
        result.update(projects=records, figure_metadata=figures,
                      manifests=check_project_manifests(errors, PROJECTS[7:]),
                      markdown=check_markdown_links(included, errors))
    except (OSError, ValueError, KeyError, SyntaxError, subprocess.CalledProcessError) as exc:
        error(errors, "validator_cannot_complete", detail=str(exc))
    result.update(status="passed" if not errors else "failed", errors=errors)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True,
                        help="Repository-relative JSON report FILE; archived records cannot be replaced")
    args = parser.parse_args()
    result = inspect()
    try:
        output = write_report(args.out, result)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    print(json.dumps({"status": result["status"], "report": relative(output),
                      "errors": result["errors"]}, ensure_ascii=False, indent=2))
    return 0 if result["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
