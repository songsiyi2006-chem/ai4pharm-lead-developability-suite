"""Read-only current-layout checks for Tasks 5-7 and optional local reproductions.

The 2026-09-20 delivery under docs/archive is historical evidence and is never
rewritten. Use validate_repository_layout.py for all seven projects and migration
byte preservation. --out names a FILE, not the historical validation directory.
"""
from pathlib import Path
import argparse
import ast
import json
import subprocess

from PIL import Image

from validate_repository_layout import (
    ROOT, PROJECTS, MANIFESTS, bounded_path, check_markdown_links,
    check_project_manifests, error, inventory, load_json, relative, sha, write_report,
)

TASKS = PROJECTS[4:]


def inspect(reproductions=None):
    errors = []
    result = {"validator": "validate_task5_7_delivery", "layout": "projects/task05-task07",
              "scope": "Current structural/provenance checks; archived delivery records remain frozen"}
    try:
        _, included = inventory()
        parsed, figures = [], []
        for project in TASKS:
            folder = ROOT / "projects" / project
            names = [name for name in sorted(included) if name.startswith(relative(folder) + "/")]
            for name in names:
                if name.endswith(".py"):
                    ast.parse((ROOT / name).read_text(encoding="utf-8-sig"), filename=name)
                    parsed.append(name)
            images = [ROOT / name for name in names if Path(name).suffix.lower() == ".png" and Path(name).name.startswith(("fig1_", "fig2_", "fig3_", "fig4_"))]
            if len(images) != 4:
                error(errors, "requested_figure_count", project=project, expected=4, actual=len(images))
            for path in images:
                try:
                    with Image.open(path) as image:
                        dpi = image.info.get("dpi", ())
                        dimensions = list(image.size)
                        if len(dpi) != 2 or not all(abs(d - 300) < 0.02 for d in dpi):
                            error(errors, "figure_dpi", path=relative(path), actual=list(dpi))
                        image.verify()
                    figures.append({"path": relative(path), "sha256": sha(path), "dpi": list(dpi), "pixels": dimensions})
                except (OSError, ValueError) as exc:
                    error(errors, "invalid_figure", path=relative(path), detail=str(exc))
        result["python_sources_parsed"] = parsed
        result["figure_metadata"] = figures
        result["markdown"] = check_markdown_links(included, errors)
        result["manifests"] = check_project_manifests(errors, TASKS)
        comparisons = []
        if reproductions:
            for project, spec, output in zip(TASKS, MANIFESTS[4:], reproductions):
                manifest_name, key, artifact_base, _ = spec
                folder = ROOT / "projects" / project
                generated = bounded_path(ROOT, output)
                if not generated.is_dir():
                    error(errors, "reproduction_directory_missing", project=project, path=relative(generated))
                    continue
                manifest = load_json(folder / manifest_name)
                for name in manifest[key]:
                    if name == "run_log.json":
                        continue
                    current = bounded_path(folder / artifact_base, name)
                    reproduced = bounded_path(generated, name)
                    if not reproduced.is_file():
                        error(errors, "reproduction_artifact_missing", project=project, artifact=name)
                        continue
                    mode = "byte_identical"
                    if project == "task05_cyp_ddi" and name == "results_summary.json":
                        old, new = load_json(current), load_json(reproduced)
                        excluded = ("script_sha256", "task4_sha256")
                        equal = {k: v for k, v in old.items() if k not in excluded} == {k: v for k, v in new.items() if k not in excluded}
                        mode = "all_fields_equal_except_historical_source_hashes"
                    else:
                        equal = sha(current) == sha(reproduced)
                    comparisons.append({"project": project, "artifact": name, "comparison": mode, "passed": equal})
                    if not equal:
                        error(errors, "independent_reproduction_mismatch", project=project, artifact=name, comparison=mode)
        result["independent_reproduction_comparisons"] = comparisons
        result["reproduction_exclusions"] = [
            "run_log.json: wall-clock timing differs",
            "manifest.json: metadata references run timing and current sources",
            "Task 5 results_summary.json: only script_sha256/task4_sha256 excluded; all other fields compared",
        ]
    except (OSError, ValueError, KeyError, SyntaxError, subprocess.CalledProcessError) as exc:
        error(errors, "validator_cannot_complete", detail=str(exc))
    result.update(status="passed" if not errors else "failed", errors=errors)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True, help="Repository-relative JSON report FILE")
    parser.add_argument("--reproductions", nargs=3, type=Path, metavar=("TASK5", "TASK6", "TASK7"))
    args = parser.parse_args()
    result = inspect(args.reproductions)
    try:
        output = write_report(args.out, result)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    print(json.dumps({"status": result["status"], "report": relative(output), "errors": result["errors"]}, ensure_ascii=False, indent=2))
    return 0 if result["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
