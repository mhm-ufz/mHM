#!/usr/bin/env python3
"""Stage coverage artifacts and verify source paths for GitLab and Coveralls."""

import argparse
import shutil
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

from lcov_cobertura.lcov_cobertura import LcovCobertura


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("build_dir", type=Path)
    parser.add_argument("--output-dir", type=Path, default=Path("coverage"))
    args = parser.parse_args()
    source_dir = Path(__file__).resolve().parents[1]
    tracked = set(
        subprocess.check_output(["git", "-C", str(source_dir), "ls-files", "-z"])
        .decode()
        .split("\0")
    )

    # Preserve filename case and original Fortran line numbers. Absolute paths
    # are needed by the converter; repository-relative LCOV is portable.
    absolute_lcov = []
    relative_lcov = []
    for line in (args.build_dir / "mhm_coverage_CI.info").read_text().splitlines():
        if line.startswith("SF:"):
            source = Path(line[3:])
            if not source.is_absolute():
                source = source_dir / source
            try:
                relative = source.resolve().relative_to(source_dir).as_posix()
            except ValueError:
                parser.error(f"Coverage source is outside the repository: {source}")
            if relative not in tracked:
                parser.error(f"Coverage source is not tracked by Git: {relative}")
            absolute_lcov.append(f"SF:{source}")
            relative_lcov.append(f"SF:{relative}")
        else:
            absolute_lcov.append(line)
            relative_lcov.append(line)

    xml = LcovCobertura("\n".join(absolute_lcov) + "\n", str(source_dir)).convert()
    report = ET.fromstring(xml)
    classes = report.findall("./packages/package/classes/class")
    if not classes or int(report.attrib["lines-valid"]) == 0:
        parser.error("Coverage report contains no executable source lines")
    for entry in classes:
        if entry.attrib["filename"] not in tracked:
            parser.error(f"Unmatched Cobertura filename: {entry.attrib['filename']}")

    shutil.copytree(
        args.build_dir / "mhm_coverage_CI", args.output_dir, dirs_exist_ok=True
    )
    (args.output_dir / "coverage.info").write_text("\n".join(relative_lcov) + "\n")
    (args.output_dir / "coverage.xml").write_text(xml)
    print(f"Coverage artifacts: {args.output_dir} ({len(classes)} source files)")
    print(f"lines......: {100 * float(report.attrib['line-rate']):.2f}%")


if __name__ == "__main__":
    main()
