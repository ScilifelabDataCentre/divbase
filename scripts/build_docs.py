"""
Auto generate documentation for DivBase CLI's typer commands.

This script is used as a hook for mkdocs (see mkdocs.yml).
In this way it ensures that the documentation is always up to date.

The script relies on typer's CLI utility to auto generate markdown documentation
https://typer.tiangolo.com/tutorial/typer-command/#generate-docs
"""

import subprocess
from pathlib import Path

SAVE_DIR = Path("docs/cli/_auto_generated")
CLI_ENTRY_POINT = Path("packages/divbase-cli/src/divbase_cli/divbase_cli.py")
CLI_COMMANDS_SRC_DIR = Path("packages/divbase-cli/src/divbase_cli/cli_commands")

SUB_COMMANDS = {
    "auth": CLI_COMMANDS_SRC_DIR / "auth_cli.py",
    "config": CLI_COMMANDS_SRC_DIR / "user_config_cli.py",
    "dimensions": CLI_COMMANDS_SRC_DIR / "dimensions_cli.py",
    "files": CLI_COMMANDS_SRC_DIR / "file_cli.py",
    "query": CLI_COMMANDS_SRC_DIR / "query_cli.py",
    "task-history": CLI_COMMANDS_SRC_DIR / "task_history_cli.py",
    "version": CLI_COMMANDS_SRC_DIR / "version_cli.py",
}


def on_startup(command, dirty):
    """
    Build CLI documentation using typer's built-in doc generation.

    (The function has to be called this to be picked up by mkdocs as a hook.)

    NOTE: An alternative hook to use could be: 'on_pre_build',
    but that triggers on every save when working on the docs which is a bit excessive.
    """
    result = subprocess.run(
        [
            "typer",
            str(CLI_ENTRY_POINT),
            "utils",
            "docs",
            "--name",
            "divbase-cli",
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    # Keep only the top level command output (e.g. divbase-cli --help)
    overview_content = result.stdout.split("##")[0].strip()
    output_file = SAVE_DIR / "divbase-cli-overview.md"
    output_file.write_text(overview_content)

    for cmd_name, cmd_path in SUB_COMMANDS.items():
        output_file = SAVE_DIR / f"{cmd_name}.md"
        subprocess.run(
            [
                "typer",
                str(cmd_path),
                "utils",
                "docs",
                "--name",
                f"divbase-cli {cmd_name}",
                "--output",
                str(output_file),
            ],
            check=True,
        )
