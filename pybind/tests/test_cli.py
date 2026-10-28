"""Installed command paths and forward-run behavior against a small routing run."""

import shutil
import subprocess
import sys
import sysconfig
from pathlib import Path

import numpy as np
import pytest
from netCDF4 import Dataset

import mhm
from mhm import cli


def console(name):
    path = shutil.which(name, path=sysconfig.get_path("scripts"))
    assert path is not None, f"missing installed console script: {name}"
    return [path]


@pytest.mark.parametrize(
    "command",
    [
        [sys.executable, "-m", "mhm"],
        console("mhm"),
        [sys.executable, "-m", "mhm.download"],
        console("mhm-download"),
    ],
)
def test_help_and_version(command):
    help_result = subprocess.run(
        command + ["--help"], capture_output=True, text=True, check=True
    )
    assert help_result.stdout.startswith("usage: mhm")
    assert not help_result.stderr
    version = subprocess.run(
        command + ["--version"], capture_output=True, text=True, check=True
    )
    expected = mhm.__version__ if "download" in " ".join(command) else mhm.f_version
    assert version.stdout.strip() == expected
    assert not version.stderr


def test_cli_arguments(monkeypatch, tmp_path):
    calls = []
    monkeypatch.setattr(cli, "run_forward", lambda **kwargs: calls.append(kwargs))
    cli.mhm(
        [
            "-n",
            "custom.nml",
            "-p",
            "params.nml",
            "-o",
            "out.nml",
            "-qq",
            "-v",
            "--log-date",
            "--log-time",
            "--log-line",
            "--log-scope",
            "#",
            "--log-file",
            "run.log",
            str(tmp_path),
        ]
    )
    assert calls == [
        dict(
            main_file="custom.nml",
            para_file="params.nml",
            out_file="out.nml",
            cwd=str(tmp_path),
            verbose=1,
            quiet=2,
            log_date=True,
            log_time=True,
            log_line=True,
            log_scope="#",
            log_file="run.log",
        )
    ]


def test_python_distribution_version():
    result = subprocess.run(
        [sys.executable, "-m", "mhm", "--pyversion"],
        capture_output=True,
        text=True,
        check=True,
    )
    assert result.stdout.strip() == mhm.__version__


def test_unknown_option_fails():
    result = subprocess.run(
        [sys.executable, "-m", "mhm", "--optimize"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 2
    assert "unrecognized arguments: --optimize" in result.stderr


@pytest.fixture
def run_root(tmp_path):
    source = Path(__file__).resolve().parents[2]
    run = tmp_path / "run directory"
    inputs = run / "examples/fixtures"
    inputs.mkdir(parents=True)
    for name in ["input.nc", "runoff.nc"]:
        shutil.copyfile(source / "examples/fixtures" / name, inputs / name)
    config = (source / "test_nml/mrm_minimal.nml").read_text()
    config = config.replace('"1990-06"', '"1990-01-03"')
    (run / "mhm.nml").write_text(config)
    shutil.copyfile(
        source / "test_nml/mhm_parameter_v6_routing.nml", run / "parameters.nml"
    )
    shutil.copyfile(source / "mhm-output-template.nml", run / "outputs.nml")
    return run


def discharge(run):
    with Dataset(run / "mrm_output.nc") as nc:
        values = nc["discharge"][:]
    assert values.count() > 0
    assert np.all(np.isfinite(values.compressed()))
    assert np.any(values.compressed() > 0.0)
    return values


@pytest.mark.parametrize("command", [[sys.executable, "-m", "mhm"], console("mhm")])
@pytest.mark.parametrize("explicit_paths", [False, True])
def test_forward_commands(run_root, command, explicit_paths, tmp_path):
    before = Path.cwd()
    mhm.run_forward(cwd=run_root, quiet=3)
    reference = discharge(run_root)
    args = ["-qqq", str(run_root)]
    if explicit_paths:
        for original, renamed in [
            ("mhm.nml", "configuration.nml"),
            ("parameters.nml", "routing.nml"),
            ("outputs.nml", "requested.nml"),
        ]:
            (run_root / original).rename(run_root / renamed)
        args += ["-n", "configuration.nml", "-p", "routing.nml", "-o", "requested.nml"]
    log = tmp_path / "forward.log"
    result = subprocess.run(
        command + args + ["--log-file", str(log)], capture_output=True, text=True
    )
    assert result.returncode == 0, result.stdout + result.stderr + log.read_text()
    assert log.is_file()
    np.testing.assert_array_equal(discharge(run_root), reference)
    assert Path.cwd() == before


def test_forward_defaults_and_repeated_logging(run_root, monkeypatch, capfd):
    monkeypatch.chdir(run_root)
    mhm.run_forward(quiet=3)
    reference = discharge(run_root)
    capfd.readouterr()
    # File units are closed/restored, and quiet settings do not accumulate across runs.
    mhm.run_forward(
        quiet=0,
        log_file=run_root / "run.log",
        log_date=True,
        log_time=True,
        log_line=True,
    )
    log = (run_root / "run.log").read_text()
    # The log prefix and complete banner row must remain on one line on every compiler.
    assert any(
        "mo_mhm_runner.F90" in line
        and "MULTISCALE HYDROLOGIC MODEL" in line
        and line.endswith("|")
        for line in log.splitlines()
    )
    mhm.run_forward(quiet=0)
    output = capfd.readouterr().out
    assert "MULTISCALE HYDROLOGIC MODEL" in output
    assert "DESTROY COMPONENTS" in output
    assert "mo_mhm_runner.F90" not in output
    np.testing.assert_array_equal(discharge(run_root), reference)
    assert Path.cwd() == run_root


def test_forward_missing_namelist_fails(tmp_path):
    result = subprocess.run(
        [sys.executable, "-m", "mhm", "-n", "missing.nml", str(tmp_path)],
        capture_output=True,
        text=True,
    )
    assert result.returncode != 0
    assert "Error reading project config:" in result.stdout + result.stderr
