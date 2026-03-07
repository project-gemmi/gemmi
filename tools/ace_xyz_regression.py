#!/usr/bin/env python3

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def load_gemmi(repo_root: Path):
    sys.path.insert(0, str(repo_root / "build" / "py"))
    import gemmi  # type: ignore
    return gemmi


def parse_args():
    repo_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(
        description="Generate ChemComp XYZ from restraints and validate with `gemmi validate -m`.")
    parser.add_argument(
        "--components",
        nargs="*",
        default=None,
        help="Component codes to test. If omitted, scan all CIFs under --ccd-dir.")
    parser.add_argument(
        "--ccd-dir",
        default=str(repo_root / "tests" / "ccd"),
        help="Directory with component CIF files.")
    parser.add_argument(
        "--recursive",
        action="store_true",
        help="Recurse under --ccd-dir when discovering CIF files.")
    parser.add_argument(
        "--limit",
        type=int,
        default=0,
        help="Limit the number of discovered components (0 means no limit).")
    parser.add_argument(
        "--gemmi-bin",
        default=str(repo_root / "build" / "gemmi"),
        help="Path to built gemmi executable.")
    parser.add_argument(
        "--mode",
        choices=["raw", "prepare"],
        default="prepare",
        help="`prepare`: run prepare_chemcomp first (recommended); `raw`: use CIF restraints as-is.")
    parser.add_argument(
        "--tables-dir",
        default=os.environ.get("GEMMI_ACEDRG_TABLES") or (
            os.path.join(os.environ["CCP4"], "share", "acedrg", "tables")
            if "CCP4" in os.environ else None),
        help="AceDRG tables directory for --mode=prepare.")
    parser.add_argument(
        "--z-score",
        type=float,
        default=2.0,
        help="Z-score threshold passed to `gemmi validate -m`.")
    parser.add_argument(
        "--keep-temp",
        action="store_true",
        help="Keep generated CIFs in a temporary directory.")
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print per-component details.")
    return parser.parse_args(), repo_root


def ensure_exists(path: Path, what: str):
    if not path.exists():
        raise SystemExit(f"{what} not found: {path}")


def build_tables(gemmi, mode: str, tables_dir: str | None):
    tables = gemmi.AcedrgTables()
    if mode == "prepare":
        if not tables_dir:
            raise SystemExit("--mode=prepare requires --tables-dir or CCP4/GEMMI_ACEDRG_TABLES")
        tables.load_tables(tables_dir)
    return tables


def generated_cif_path(tmpdir: Path, comp_id: str) -> Path:
    return tmpdir / f"{comp_id}_generated.cif"


def resolve_component_paths(ccd_dir: Path, components, recursive: bool):
    if components:
        paths = []
        for comp_id in components:
            path = ccd_dir / f"{comp_id}.cif"
            ensure_exists(path, "component CIF")
            paths.append(path)
        return paths
    pattern = "**/*.cif" if recursive else "*.cif"
    return sorted(ccd_dir.glob(pattern))


def run_component(gemmi, gemmi_bin: Path, source: Path, mode: str,
                  tables, z_score: float, out_dir: Path):
    comp_id = source.stem

    doc = gemmi.cif.read(str(source))
    cc = gemmi.make_chemcomp_from_block(doc.sole_block())
    if mode == "prepare":
        gemmi.prepare_chemcomp(cc, tables)

    for atom in cc.atoms:
        atom.xyz = gemmi.Position(float("nan"), float("nan"), float("nan"))

    placed = gemmi.generate_chemcomp_xyz_from_restraints(cc)
    generated = generated_cif_path(out_dir, comp_id)
    out_doc = gemmi.cif.Document()
    out_block = out_doc.add_new_block("comp_" + cc.name)
    gemmi.add_chemcomp_to_block(cc, out_block, [], False)
    generated.write_text(out_doc.as_string(), encoding="utf-8")

    cmd = [str(gemmi_bin), "validate", "-m", f"--z-score={z_score}", str(generated)]
    proc = subprocess.run(cmd, text=True, capture_output=True, check=False)
    output = (proc.stdout + proc.stderr).strip()
    ok = proc.returncode == 0 and output == "" and placed == len(cc.atoms)
    return {
        "component": comp_id,
        "placed": placed,
        "atoms": len(cc.atoms),
        "returncode": proc.returncode,
        "output": output,
        "generated": str(generated),
        "ok": ok,
    }


def main():
    args, repo_root = parse_args()
    gemmi = load_gemmi(repo_root)
    ccd_dir = Path(args.ccd_dir)
    gemmi_bin = Path(args.gemmi_bin)
    ensure_exists(ccd_dir, "CCD directory")
    ensure_exists(gemmi_bin, "gemmi executable")
    tables = build_tables(gemmi, args.mode, args.tables_dir)
    sources = resolve_component_paths(ccd_dir, args.components, args.recursive)
    if args.limit > 0:
        sources = sources[:args.limit]
    if not sources:
        raise SystemExit(f"no CIF files found under {ccd_dir}")

    out_dir = Path(tempfile.mkdtemp(prefix="ace_xyz_"))
    failures = []
    results = []
    try:
        for source in sources:
            result = run_component(gemmi, gemmi_bin, source, args.mode,
                                   tables, args.z_score, out_dir)
            results.append(result)
            if args.verbose or not result["ok"]:
                status = "OK" if result["ok"] else "FAIL"
                print(f"{status} {result['component']}: placed {result['placed']}/{result['atoms']}")
                if result["output"]:
                    print(result["output"])
                if not result["ok"]:
                    print(f"generated: {result['generated']}")
            if not result["ok"]:
                failures.append(result["component"])
    finally:
        if args.keep_temp:
            print(f"kept generated CIFs in {out_dir}")
        else:
            shutil.rmtree(out_dir, ignore_errors=True)

    if failures:
        print(f"FAILED components: {', '.join(failures)}", file=sys.stderr)
        return 1

    print(f"All {len(results)} components passed in {args.mode} mode.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
