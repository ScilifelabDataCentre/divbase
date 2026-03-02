"""
Module responsible for validating the CLI version from incoming requests to the API.

Can be used to outright reject a request (via Middleware) from a CLI version that is too old
OR
On login, add annoucnement that a new version of the CLI can be installed.
"""

from divbase_api.api_config import settings


def _parse_version(version: str) -> list[int]:
    """
    Parse a semantic version string like '0.1.0' '2.1.0-rc123' or '3.1.0.dev3'
    into [major, minor, patch].
    Things like dev20, -alpha5, beta3, rc12 etc.. are ignored.
    """
    parts = version.split(".")[:3]
    numeric_parts = []
    for idx, part in enumerate(parts):
        if idx < 2:
            # major and minor: always just convert to int
            numeric_parts.append(int(part))
        else:
            # patch: take leading digits only
            num = ""
            for c in part:
                if c.isdigit():
                    num += c
                else:
                    break
            numeric_parts.append(int(num) if num else 0)
    return numeric_parts


def cli_version_outdated(cli_version: str, minimum_cli_version: str = settings.api.minimum_cli_version) -> bool:
    """
    Compare the CLI version against the minimum required version.
    Version strings are expected to be in the format "major.minor.patch" (e.g., "2.5.0").
    Any extra parts, e.g. alpha, beta, rc1, dev3 etc.. ignored.
    """
    cli_parts = _parse_version(cli_version)
    min_parts = _parse_version(minimum_cli_version)

    for cli_part, min_part in zip(cli_parts, min_parts, strict=True):
        if cli_part > min_part:
            return False
        elif cli_part < min_part:
            return True

    return False  # versions are equal


def cli_update_available(cli_version: str, latest_cli_version: str = settings.api.latest_cli_version) -> bool:
    """
    Compare the CLI version against the latest available version to determine if an update is available.
    Version strings are expected to be in the format "major.minor.patch" (e.g., "2.5.0").
    Any extra parts, e.g. alpha, beta, rc1, dev3 etc.. ignored.
    """
    cli_parts = _parse_version(cli_version)
    latest_parts = _parse_version(latest_cli_version)

    for cli_part, latest_part in zip(cli_parts, latest_parts, strict=True):  # noqa: SIM110, (more readable)
        if cli_part < latest_part:
            return True
        elif cli_part > latest_part:
            return False  # early exit as e.g. 1.0.0 > 0.9.0

    return False  # versions are equal or cli version is newer than server knows about
