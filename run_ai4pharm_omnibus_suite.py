#!/usr/bin/env python3
"""Ten-task local CPU runner, traceable figure rebuild and offline navigation.

Archived mode reuses published calculations. Recompute mode executes all ten
original drivers sequentially in a NEW output tree. No automatic clinical
parameter transfer, network downloads, environment installation or Git writes.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import signal
import subprocess
import sys
import time

# Set before importing any numerical package, including figure helpers.
for _key in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_key] = "1"
os.environ["MPLBACKEND"] = "Agg"

from omnibus.catalog import PROJECTS, DRIVERS, TITLES, FIGURES, BOUNDARIES

ROOT = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False)+"\n", encoding="utf-8", newline="\n")


def prepare_output(out, repo=ROOT):
    """Refuse repository/input overwrite, including resolved symlink aliases."""
    out, repo = Path(out).resolve(), Path(repo).resolve()
    if out == repo or repo.is_relative_to(out):
        raise ValueError("Output must not be the repository or its ancestor")
    protected = ("projects", "docs", "omnibus", "tests", "tools", ".git", "figures_omnibus", "data_omnibus")
    if any(out.is_relative_to(repo/p) for p in protected):
        raise ValueError("Output overlaps protected repository artifacts; use work/NEW_RUN")
    if out.exists():
        raise FileExistsError("Output already exists; select a NEW directory")
    out.mkdir(parents=True)
    return out


def child_command(task, out, python=sys.executable, rna_library=None):
    if not 1 <= task <= 10:
        raise ValueError("task must be 1..10")
    target = Path(out)/"projects"/PROJECTS[task-1]
    if task in (5, 10):
        target = target/"results"
    flag = "--output-dir" if task <= 3 else "--out"
    command = [str(python), str(ROOT/"projects"/PROJECTS[task-1]/DRIVERS[task-1]), flag, str(target)]
    if task == 1:
        command += ["--conformer-timeout", "120"]
    if task == 10 and rna_library:
        command += ["--rna-library", str(Path(rna_library).resolve())]
    return command, target


def stop_process_tree(process):
    """Terminate this driver's descendants too (Task 1 has geometry workers)."""
    if os.name == "nt":
        killed = subprocess.run(["taskkill", "/PID", str(process.pid), "/T", "/F"],
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=False)
        detail = {"method": "taskkill /T /F", "returncode": killed.returncode}
        if process.poll() is None:
            process.kill()
    else:
        try:
            os.killpg(process.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        detail = {"method": "SIGKILL process group"}
    process.wait()
    return detail


def run_child(command, log_handle, timeout):
    options = {"creationflags": subprocess.CREATE_NEW_PROCESS_GROUP} if os.name == "nt" else {"start_new_session": True}
    with subprocess.Popen(command, cwd=ROOT, stdout=log_handle, stderr=subprocess.STDOUT,
                          env=os.environ.copy(), **options) as process:
        try:
            return process.wait(timeout=timeout)
        except subprocess.TimeoutExpired as exc:
            exc.tree_termination = stop_process_tree(process)
            raise
        except KeyboardInterrupt:
            stop_process_tree(process)
            raise


def recompute(out, timeout, rna_library=None):
    source = out/"calculations"
    records = []
    logs = out/"logs"
    logs.mkdir()
    for task in range(1, 11):
        command, target = child_command(task, source, rna_library=rna_library)
        target.parent.mkdir(parents=True, exist_ok=True)
        started = time.perf_counter()
        print(f"[{task:02d}/10] Recomputing {PROJECTS[task-1]} (one numerical thread)", flush=True)
        log = logs/f"task{task:02d}.txt"
        status, returncode, termination = "running", None, None
        try:
            with log.open("w", encoding="utf-8", newline="\n") as handle:
                returncode = run_child(command, handle, timeout)
            status = "passed" if returncode == 0 else "failed"
        except subprocess.TimeoutExpired as exc:
            status = "timeout"
            termination = getattr(exc, "tree_termination", None)
        record = {"task": task, "status": status, "returncode": returncode,
                  "elapsed_seconds": time.perf_counter()-started,
                  "command": command, "log": log.relative_to(out).as_posix(),
                  "process_tree_termination": termination,
                  "driver_sha256": sha(ROOT/"projects"/PROJECTS[task-1]/DRIVERS[task-1])}
        records.append(record)
        write_json(out/"data_omnibus"/"execution.json", records)
        if status != "passed":
            raise RuntimeError(f"Task {task} {status}; preserved partial outputs and log: {log}")
        print(f"[{task:02d}/10] Completed in {record['elapsed_seconds']:.1f} s", flush=True)
    return source, records


def verify_archives():
    sys.path.insert(0, str(ROOT/"tools"))
    from validate_repository_layout import check_project_manifests
    errors = []
    records = check_project_manifests(errors)
    if errors:
        raise ValueError(f"Archive integrity check failed: {errors[:3]}")
    return records


def versions():
    result = {"python": platform.python_version(), "platform": platform.platform()}
    for name in ("numpy", "scipy", "matplotlib", "rdkit", "pandas", "ViennaRNA"):
        try:
            result[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            result[name] = "not_installed_in_main_environment"
    return result


def rebuild(source, out, mode, execution):
    from omnibus.figures_1_5 import generate as first
    from omnibus.figures_6_10 import generate as second
    from omnibus.presentation import write_navigation
    from PIL import Image

    figures = out/"figures_omnibus"
    figures.mkdir()
    (out/"data_omnibus").mkdir(exist_ok=True)
    records = first(source, figures)+second(source, figures)
    if sorted(r["task"] for r in records) != list(range(1, 11)):
        raise ValueError("Figure generators must return exactly the ten distinct tasks")
    inputs = {}
    for record in records:
        task = record["task"]
        if record["file"] != FIGURES[task-1]:
            raise ValueError(f"Unexpected figure filename for Task {task}")
        path = figures/record["file"]
        with Image.open(path) as im:
            dpi = im.info.get("dpi", ())
            if len(dpi) != 2 or any(abs(d-300) > .02 for d in dpi):
                raise ValueError(f"Not 300 DPI: {path}")
            record["pixels"] = list(im.size)
            record["dpi"] = list(dpi)
            im.verify()
        record["sha256"] = sha(path)
        record["evidence_boundary"] = BOUNDARIES[task-1]
        for name in record["sources"]:
            path = (source/name).resolve()
            if not path.is_relative_to(source.resolve()) or not path.is_file():
                raise ValueError(f"Invalid figure source: {name}")
            inputs[name] = sha(path)
    write_json(out/"data_omnibus"/"figure_provenance.json", records)
    write_json(out/"data_omnibus"/"source_hashes.json", inputs)
    write_navigation(out, records, ROOT, mode, source)
    info = {"status": "completed", "mode": mode,
            "semantics": "Ten separately parameterized studies; only Task 4 to 5 is an executed inter-task model link",
            "source_root": "repository_root" if source.resolve() == ROOT else str(source.resolve()), "environment": versions(),
            "numerical_threads": 1, "execution": execution, "figure_count": len(records),
            "unique_figure_input_files": len(inputs), "clinical_validation": False}
    write_json(out/"data_omnibus"/"run_summary.json", info)
    sources = [Path(__file__), *sorted((ROOT/"omnibus").glob("*.py"))]
    manifest = {"schema": "ai4pharm-omnibus-1", "scientific_driver_sha256":
                {f"projects/{p}/{d}": sha(ROOT/"projects"/p/d) for p, d in zip(PROJECTS, DRIVERS)},
                "source_code_sha256":
                {p.relative_to(ROOT).as_posix(): sha(p) for p in sources},
                "files_sha256": {p.relative_to(out).as_posix(): sha(p) for p in sorted(out.rglob("*"))
                                  if p.is_file() and "calculations" not in p.relative_to(out).parts
                                  and "__pycache__" not in p.parts}}
    write_json(out/"manifest_omnibus.json", manifest)
    return info


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True, help="NEW output directory; use work/NAME")
    parser.add_argument("--mode", choices=("archived", "recompute"), default="archived")
    parser.add_argument("--task-timeout", type=int, default=1800, help="Seconds per recomputed task")
    parser.add_argument("--rna-library", type=Path, help="Existing optional ViennaRNA library for Task 10")
    parser.add_argument("--calculations-only", action="store_true", help="Recompute and save execution records without rendering")
    parser.add_argument("--source-root", type=Path, help="Existing canonical projects tree from a previous recompute; no archived fallback")
    args = parser.parse_args(argv)
    if args.task_timeout <= 0:
        parser.error("--task-timeout must be positive")
    if args.calculations_only and args.mode != "recompute":
        parser.error("--calculations-only requires --mode recompute")
    if args.source_root and args.mode != "archived":
        parser.error("--source-root cannot be combined with recompute")
    if args.source_root and not (args.source_root/"projects").is_dir():
        parser.error("--source-root must contain projects/")
    verify_archives()
    out = prepare_output(args.out)
    execution = []
    source = args.source_root.resolve() if args.source_root else ROOT
    mode = "provided_source_tree" if args.source_root else args.mode
    try:
        if args.mode == "recompute":
            source, execution = recompute(out, args.task_timeout, args.rna_library)
        if args.calculations_only:
            write_json(out/"data_omnibus"/"run_summary.json", {
                "status": "calculations_completed", "mode": args.mode,
                "execution": execution, "environment": versions(), "figures_generated": False})
        else:
            result = rebuild(source, out, mode, execution)
            print(json.dumps({k: result[k] for k in ("status", "mode", "figure_count", "unique_figure_input_files")}), flush=True)
    except Exception as exc:
        write_json(out/"data_omnibus"/"failure.json", {"status": "failed", "type": type(exc).__name__, "detail": str(exc)})
        raise
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
