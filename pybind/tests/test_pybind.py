"""Tests for the temporary version-only package and example downloader."""

import io
import tarfile
from urllib.parse import parse_qs, unquote, urlparse

import pytest

import mhm
from mhm import download


def test_versions():
    assert isinstance(mhm.__version__, str)
    assert isinstance(mhm.f_version, str)
    assert mhm.f_version.strip() == mhm.f_version
    assert mhm.f_version
    for removed in ["model", "run", "get", "set", "set_meteo", "get_variable"]:
        with pytest.raises(AttributeError):
            getattr(mhm, removed)


@pytest.mark.parametrize("domain", [1, 2])
@pytest.mark.parametrize(
    "version,branch", [("6.0.0.dev1", "v6alpha"), ("6.0.0", "v6.0.0")]
)
def test_download(tmp_path, monkeypatch, domain, version, branch):
    monkeypatch.setattr(mhm, "__version__", version)
    name = f"domain_{domain:02d}"

    def retrieve(url, archive):
        assert unquote(urlparse(url).path).endswith(f"/{branch}/mhm-{branch}.tar.gz")
        assert parse_qs(urlparse(url).query)["path"] == [f"examples/{name}"]
        with tarfile.open(archive, "w:gz") as tar:
            for leaf in [
                "mhm.nml",
                "parameters.nml",
                "outputs.nml",
                "input/meteo/pre.nc",
            ]:
                entry = tarfile.TarInfo(f"mhm-{branch}/examples/{name}/{leaf}")
                entry.size = 7
                tar.addfile(entry, io.BytesIO(b"fixture"))

    monkeypatch.setattr(download, "urlretrieve", retrieve)
    monkeypatch.chdir(tmp_path)
    mhm.download_examples(domain=domain)
    assert (tmp_path / name / "input/meteo/pre.nc").read_bytes() == b"fixture"
    assert (tmp_path / name / "parameters.nml").is_file()


def test_explicit_branch_path_and_cli(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(download, "download_examples", lambda *args: calls.append(args))
    download.cli(
        ["--branch", "my/ref", "--domain", "2", "--path", str(tmp_path / "custom")]
    )
    assert calls == [("my/ref", 2, str(tmp_path / "custom"), False)]


def test_invalid_domain_does_not_download(monkeypatch):
    monkeypatch.setattr(
        download, "urlretrieve", lambda *args: pytest.fail("unexpected download")
    )
    with pytest.raises(ValueError, match="needs to be 1 or 2"):
        mhm.download_examples(domain=3)


@pytest.mark.parametrize(
    "entry", ["mhm-ref/examples/domain_01/../escape", "mhm-ref/unrelated.txt"]
)
def test_invalid_archive(tmp_path, monkeypatch, entry):
    def retrieve(url, archive):
        with tarfile.open(archive, "w:gz") as tar:
            item = tarfile.TarInfo(entry)
            item.size = 1
            tar.addfile(item, io.BytesIO(b"x"))

    monkeypatch.setattr(download, "urlretrieve", retrieve)
    with pytest.raises(ValueError):
        mhm.download_examples(path=tmp_path / "target")
    assert not (tmp_path / "target").exists()
