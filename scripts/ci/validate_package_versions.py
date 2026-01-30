"""
Used in CI to validate that package versions are correctly set before publishing to PyPI.

Catches 2 cases:

1. We keep all 3 packages in this monorepo at the same version for simplicity, so all packages should have the same version in their __init__.py files.
   This is to ensure that the divbase-cli, divbase-lib, and divbase-api packages all have the same version number set.

2. The divbase-cli package depends on divbase-lib and its dependency version should match the actual version of divbase-lib.
    So the version of divbase-lib defined in divbase-cli/pyproject.toml should match the version defined in divbase-lib/__init__.py.

NOTE: divbase-api is not published as a package so we don't need to fix the divbase-lib version in the pyproject.toml,
instead it just uses the one in the repo when it the image is built.

NOTE: As we don't want to install the packages in the CI environment just to check versions,
we read the versions from the source files directly.
"""

import re
from pathlib import Path

import tomllib

ROOT_PATH = Path(__file__).parent.parent.parent
DIVBASE_LIB_INIT_PATH = ROOT_PATH / "packages/divbase-lib/src/divbase_lib/__init__.py"
DIVBASE_CLI_INIT_PATH = ROOT_PATH / "packages/divbase-cli/src/divbase_cli/__init__.py"
DIVBASE_API_INIT_PATH = ROOT_PATH / "packages/divbase-api/src/divbase_api/__init__.py"
DIVBASE_CLI_PYPROJECT_PATH = ROOT_PATH / "packages/divbase-cli/pyproject.toml"


def get_version_from_init(file_path: Path) -> str:
    content = file_path.read_text()
    match = re.search(r'__version__ = "([^"]+)"', content)
    if not match:
        raise ValueError(f"Could not find __version__ in {file_path}")
    return match.group(1)


if __name__ == "__main__":
    divbase_lib_version = get_version_from_init(DIVBASE_LIB_INIT_PATH)
    divbase_cli_version = get_version_from_init(DIVBASE_CLI_INIT_PATH)
    divbase_api_version = get_version_from_init(DIVBASE_API_INIT_PATH)

    if not (divbase_lib_version == divbase_cli_version == divbase_api_version):
        raise ValueError(
            f"Version mismatch detected between packages! \n"
            f"divbase-lib version: {divbase_lib_version} \n"
            f"divbase-cli version: {divbase_cli_version} \n"
            f"divbase-api version: {divbase_api_version} \n"
            "All 3 package versions should be the same."
        )

    with open(DIVBASE_CLI_PYPROJECT_PATH, "rb") as f:
        py_project_data = tomllib.load(f)
    dependencies = py_project_data.get("project", {}).get("dependencies", [])

    cli_dep_version = None
    for dep in dependencies:
        if dep.startswith("divbase-lib=="):
            cli_dep_version = dep.split("==")[1]

    if not cli_dep_version:
        raise ValueError("Could not find divbase-lib dependency in divbase-cli/pyproject.toml")

    if divbase_lib_version != cli_dep_version:
        raise ValueError(
            f"Version mismatch! divbase-lib version: {divbase_lib_version} \n"
            f"divbase-cli pyproject.toml's dependency of divbase-lib version: {cli_dep_version} \n"
            "It is not expected that these versions should differ."
        )

    print(f"Divbase package versions validated successfully, current version: {divbase_lib_version}")
