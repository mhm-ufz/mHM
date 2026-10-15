"""Compiled metadata and deterministic, offline downloader tests."""

import io
import tarfile
from urllib.error import ContentTooShortError, HTTPError, URLError
from urllib.parse import parse_qs, unquote, urlparse

import pytest

import mhm
from mhm import download


def archive_at(path, domain=1, entries=None):
    """Create a repository-shaped archive, including an unrelated file."""
    if entries is None:
        prefix = f"mhm-ref/examples/domain_{domain:02d}"
        entries = [
            "mhm-ref/README.md",
            *[
                f"{prefix}/{leaf}"
                for leaf in [
                    "mhm.nml",
                    "parameters.nml",
                    "outputs.nml",
                    "input/meteo/pre.nc",
                ]
            ],
        ]
    with tarfile.open(path, "w:gz") as tar:
        for name in entries:
            entry = tarfile.TarInfo(name)
            entry.size = 7
            tar.addfile(entry, io.BytesIO(b"fixture"))


def test_versions():
    assert isinstance(mhm.__version__, str)
    assert isinstance(mhm.f_version, str)
    assert mhm.f_version.strip() == mhm.f_version
    assert mhm.f_version
    assert callable(mhm.run_forward)
    for removed in ["model", "run", "get", "set", "set_meteo", "get_variable"]:
        with pytest.raises(AttributeError):
            getattr(mhm, removed)


@pytest.mark.parametrize("domain", [1, 2])
@pytest.mark.parametrize("source", ["auto", "github", "git.ufz"])
@pytest.mark.parametrize(
    "version,branch", [("6.0.0.dev1", "v6alpha"), ("6.0.0", "v6.0.0")]
)
def test_download(tmp_path, monkeypatch, domain, source, version, branch):
    monkeypatch.setattr(mhm, "__version__", version)
    name = f"domain_{domain:02d}"
    calls = []

    def retrieve(url, archive):
        calls.append(url)
        parsed = urlparse(url)
        if source == "git.ufz":
            assert parsed.hostname == "git.ufz.de"
            assert unquote(parsed.path).endswith(f"/{branch}/mhm-{branch}.tar.gz")
            assert parse_qs(parsed.query)["path"] == [f"examples/{name}"]
        else:
            assert parsed.hostname == "github.com"
            assert unquote(parsed.path).endswith(f"/archive/{branch}.tar.gz")
        archive_at(archive, domain)

    monkeypatch.setattr(download, "urlretrieve", retrieve)
    monkeypatch.chdir(tmp_path)
    mhm.download_examples(domain=domain, source=source)
    assert len(calls) == 1
    assert (tmp_path / name / "input/meteo/pre.nc").read_bytes() == b"fixture"
    assert (tmp_path / name / "parameters.nml").is_file()
    assert not (tmp_path / name / "README.md").exists()


def test_explicit_branch_path_and_cli(tmp_path, monkeypatch):
    calls = []
    monkeypatch.setattr(download, "download_examples", lambda *args: calls.append(args))
    download.cli(
        [
            "--branch",
            "my/ref",
            "--domain",
            "2",
            "--path",
            str(tmp_path / "custom"),
            "--source",
            "github",
        ]
    )
    assert calls == [("my/ref", 2, str(tmp_path / "custom"), False, "github")]


@pytest.mark.parametrize("source", ["github", "git.ufz"])
def test_ref_encoding(tmp_path, monkeypatch, source):
    def retrieve(url, archive):
        assert "feature%2Fnew%20examples" in url
        archive_at(archive)

    monkeypatch.setattr(download, "urlretrieve", retrieve)
    mhm.download_examples(
        branch="feature/new examples", source=source, path=tmp_path / "target"
    )


@pytest.mark.parametrize(
    "failure", ["http", "network", "partial", "tar", "missing", "unsafe"]
)
def test_github_falls_back_to_ufz(tmp_path, monkeypatch, failure):
    calls = []

    def retrieve(url, archive):
        host = urlparse(url).hostname
        calls.append(host)
        if host == "github.com":
            if failure == "http":
                raise HTTPError(url, 404, "missing ref", None, None)
            if failure == "network":
                raise URLError("unreachable")
            if failure == "partial":
                raise ContentTooShortError("incomplete download", None)
            if failure == "tar":
                archive.write_bytes(b"not an archive")
            if failure == "missing":
                archive_at(archive, domain=2)
            if failure == "unsafe":
                archive_at(
                    archive,
                    entries=[
                        "mhm-ref/examples/domain_01/stale.txt",
                        "mhm-ref/examples/domain_01/../escape",
                    ],
                )
        else:
            archive_at(archive)

    monkeypatch.setattr(download, "urlretrieve", retrieve)
    mhm.download_examples(path=tmp_path / "target")
    assert calls == ["github.com", "git.ufz.de"]
    assert (tmp_path / "target/mhm.nml").read_bytes() == b"fixture"
    assert not (tmp_path / "target/stale.txt").exists()


@pytest.mark.parametrize(
    "source,host", [("github", "github.com"), ("git.ufz", "git.ufz.de")]
)
def test_explicit_source_failure_does_not_fall_back(
    tmp_path, monkeypatch, source, host
):
    calls = []

    def retrieve(url, archive):
        calls.append(urlparse(url).hostname)
        raise URLError("unreachable")

    monkeypatch.setattr(download, "urlretrieve", retrieve)
    with pytest.raises(RuntimeError, match=source):
        mhm.download_examples(path=tmp_path / "target", source=source)
    assert calls == [host]


def test_all_sources_fail(tmp_path, monkeypatch):
    def retrieve(url, archive):
        raise URLError("unreachable")

    monkeypatch.setattr(download, "urlretrieve", retrieve)
    with pytest.raises(RuntimeError, match=r"Attempts: github:.*git.ufz:"):
        mhm.download_examples(path=tmp_path / "target")
    assert not (tmp_path / "target").exists()


@pytest.mark.parametrize("kwargs", [{"domain": 3}, {"source": "unknown"}])
def test_invalid_selection_does_not_download(monkeypatch, kwargs):
    monkeypatch.setattr(
        download, "urlretrieve", lambda *args: pytest.fail("unexpected download")
    )
    with pytest.raises(ValueError, match="needs to be"):
        mhm.download_examples(**kwargs)


@pytest.mark.parametrize(
    "entries",
    [
        ["mhm-ref/examples/domain_01/../escape"],
        ["/mhm-ref/examples/domain_01/escape"],
        ["mhm-ref/examples/domain_01/..\\escape"],
        ["mhm-ref/unrelated.txt"],
        ["root-a/examples/domain_01/a", "root-b/examples/domain_01/b"],
    ],
)
def test_invalid_archive(tmp_path, monkeypatch, entries):
    monkeypatch.setattr(
        download,
        "urlretrieve",
        lambda url, archive: archive_at(archive, entries=entries),
    )
    with pytest.raises(RuntimeError, match="failed to download"):
        mhm.download_examples(path=tmp_path / "target")
    assert not (tmp_path / "target").exists()
    assert not (tmp_path / "escape").exists()


def test_archive_link_is_not_extracted(tmp_path, monkeypatch):
    def retrieve(url, archive):
        with tarfile.open(archive, "w:gz") as tar:
            entry = tarfile.TarInfo("mhm-ref/examples/domain_01/link")
            entry.type = tarfile.SYMTYPE
            entry.linkname = "../../outside"
            tar.addfile(entry)

    monkeypatch.setattr(download, "urlretrieve", retrieve)
    with pytest.raises(RuntimeError, match="Unsupported entry"):
        mhm.download_examples(path=tmp_path / "target")
    assert not (tmp_path / "target").exists()


def test_local_destination_error_does_not_fall_back(tmp_path, monkeypatch):
    calls = []

    def retrieve(url, archive):
        calls.append(urlparse(url).hostname)
        archive_at(archive)

    monkeypatch.setattr(download, "urlretrieve", retrieve)
    target = tmp_path / "existing"
    target.mkdir()
    (target / "keep").write_text("unchanged")
    with pytest.raises(FileExistsError):
        mhm.download_examples(path=target)
    assert calls == ["github.com"]
    assert (target / "keep").read_text() == "unchanged"
