"""
Helper script to set the version of all divbase packages.
Includes setting the version of divbase-lib used by divbase-cli.

Usage:
uv run scripts/set_divbase_version.py [NEW_VERSION]
"""

import re
import sys
from pathlib import Path

ROOT_PATH = Path(__file__).parent.parent
FILES_TO_BUMP = [
    ROOT_PATH / "packages/divbase-lib/src/divbase_lib/__init__.py",
    ROOT_PATH / "packages/divbase-cli/src/divbase_cli/__init__.py",
    ROOT_PATH / "packages/divbase-api/src/divbase_api/__init__.py",
]
CLI_PYPROJECT = ROOT_PATH / "packages/divbase-cli/pyproject.toml"


def change_versions(new_version: str) -> None:
    for file_path in FILES_TO_BUMP:
        content = file_path.read_text()
        new_content = re.sub(r'__version__ = "[^"]+"', f'__version__ = "{new_version}"', content)
        file_path.write_text(new_content)


def update_cli_dependency_to_divbase_lib(new_version: str) -> None:
    content = CLI_PYPROJECT.read_text()
    new_content = re.sub(r'divbase-lib==[^" \n]+', f"divbase-lib=={new_version}", content)
    CLI_PYPROJECT.write_text(new_content)


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: uv run scripts/set_divbase_version.py [NEW_VERSION]")
        sys.exit(1)

    version_input = sys.argv[1]
    if version_input.startswith("v"):
        version_input = version_input[1:]
        print("Removed leading 'v' from your specified version input.")

    change_versions(version_input)
    update_cli_dependency_to_divbase_lib(version_input)
    print(f"divbase packages version's set to: {version_input}")
