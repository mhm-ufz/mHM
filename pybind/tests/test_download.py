"""Offline checks for legacy development and release example downloads."""

import io
import tarfile
from urllib.parse import parse_qs, unquote, urlparse

import pytest

import mhm
from mhm import download


@pytest.mark.parametrize("domain,folder", [(1, "test_domain"), (2, "test_domain_2")])
@pytest.mark.parametrize("source", ["auto", "git.ufz", "github"])
@pytest.mark.parametrize(
    "version,explicit,expected",
    [
        ("5.13.5.dev1", None, "v5.x"),
        ("5.13.4", None, "v5.13.4"),
        ("5.13.5.dev1", "feature/legacy", "feature/legacy"),
    ],
)
def test_download_ref(
    tmp_path, monkeypatch, domain, folder, source, version, explicit, expected
):
    monkeypatch.setattr(mhm, "__version__", version)
    calls = []

    def retrieve(url, archive):
        calls.append(url)
        parsed = urlparse(url)
        assert expected in unquote(parsed.path)
        if source == "github":
            assert parsed.hostname == "github.com"
        else:
            assert parsed.hostname == "git.ufz.de"
            assert parse_qs(parsed.query)["path"] == [folder]
        with tarfile.open(archive, "w:gz") as tar:
            entry = tarfile.TarInfo(f"mhm-ref/{folder}/mhm.nml")
            entry.size = 7
            tar.addfile(entry, io.BytesIO(b"fixture"))

    monkeypatch.setattr(download, "urlretrieve", retrieve)
    target = tmp_path / "domain"
    download.download_test(branch=explicit, domain=domain, path=target, source=source)
    assert len(calls) == 1
    assert (target / "mhm.nml").read_bytes() == b"fixture"
