"""
Used in CI to validate that package versions are correctly set before publishing to PyPI.

This is to catch a case where the divbase-lib version is bumped,
but the version of divbase-lib to use in divbase-cli (defined in pyproject.toml) is not also bumped.

NOTE: As we don't want to install the packages in the CI environment just to check versions,
we read the versions from the source files directly.
"""

import re
from pathlib import Path

import tomllib

ROOT_PATH = Path(__file__).parent.parent.parent
DIVBASE_LIB_INIT_PATH = ROOT_PATH / "packages/divbase-lib/src/divbase_lib/__init__.py"
DIVBASE_CLI_PYPROJECT_PATH = ROOT_PATH / "packages/divbase-cli/pyproject.toml"


if __name__ == "__main__":
    lib_init_content = DIVBASE_LIB_INIT_PATH.read_text()
    match = re.search(r'__version__ = "([^"]+)"', lib_init_content)
    if not match:
        raise ValueError("Could not find __version__ in divbase-lib/__init__.py")
    divbase_lib_version = match.group(1)

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

    print(f"Success: Versions used in divbase-lib and divbase-cli match: {divbase_lib_version}")
