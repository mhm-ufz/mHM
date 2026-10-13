"""Download runnable v6 example domains from the mHM repository."""

import argparse
import shutil
import tarfile
from pathlib import Path
from tempfile import TemporaryDirectory
from urllib.parse import quote, urlencode
from urllib.request import urlretrieve

MHM_URL = "https://git.ufz.de/mhm/mhm/-/archive/{branch}/mhm-{branch}.tar.gz?{query}"
VALID_DOMAINS = [1, 2]


def download_examples(branch=None, domain=1, path=None, verbose=False):
    """Download one example domain, preserving its runnable directory layout.

    Development packages use v6alpha; releases use their matching version tag.
    Explicit branch/tag/commit and destination arguments override these defaults.
    """
    from . import __version__

    if domain not in VALID_DOMAINS:
        raise ValueError(f"mhm-download: 'domain' needs to be 1 or 2. Got: '{domain}'")
    if branch is None:
        branch = "v6alpha" if "dev" in __version__ else f"v{__version__}"
    name = f"domain_{domain:02d}"
    folder = f"examples/{name}"
    path = Path(path or name)
    url = MHM_URL.format(
        branch=quote(branch, safe=""), query=urlencode({"path": folder})
    )
    if verbose:
        print(
            f"Downloading mHM example domain {domain} from {branch} to {path.resolve()}"
        )
    with TemporaryDirectory() as tmp_dir:
        archive = Path(tmp_dir) / "example.tar.gz"
        extracted = Path(tmp_dir) / name
        extracted.mkdir()
        urlretrieve(url, archive)
        found = False
        with tarfile.open(archive, "r:gz") as tar:
            for member in tar.getmembers():
                # GitLab archives wrap the requested repository path in a root folder.
                parts = Path(member.name).parts
                if len(parts) < 3 or parts[1:3] != ("examples", name):
                    continue
                relative = Path(*parts[3:])
                if ".." in relative.parts or relative.is_absolute():
                    raise ValueError("Invalid path in example archive")
                destination = extracted / relative
                if member.isdir():
                    destination.mkdir(parents=True, exist_ok=True)
                elif member.isfile():
                    destination.parent.mkdir(parents=True, exist_ok=True)
                    with tar.extractfile(member) as source:
                        with destination.open("wb") as target:
                            shutil.copyfileobj(source, target)
                    found = True
                else:
                    raise ValueError("Unsupported entry in example archive")
        if not found:
            raise ValueError(f"Archive does not contain {folder}")
        shutil.copytree(extracted, path)


def cli(argv=None):
    """Command line interface for example-domain downloads."""
    from . import __version__

    parser = argparse.ArgumentParser(
        description="Download a runnable mHM v6 example domain.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("-V", "--version", action="version", version=__version__)
    parser.add_argument("-v", "--verbose", action="store_true")
    parser.add_argument("-b", "--branch", help="Repository branch, tag, or commit")
    parser.add_argument("-d", "--domain", type=int, default=1, choices=VALID_DOMAINS)
    parser.add_argument(
        "-p", "--path", help="Destination; defaults to domain_01 or domain_02"
    )
    args = parser.parse_args(argv)
    download_examples(args.branch, args.domain, args.path, args.verbose)
