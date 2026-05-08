"""
Auto generate documentation for DivBase CLI's typer commands and our GitHub releases.

This script is used as a hook for mkdocs (see mkdocs.yml).
In this way it ensures that the documentation is always up to date.

The script relies on typer's CLI utility to auto generate markdown documentation
https://typer.tiangolo.com/tutorial/typer-command/#generate-docs

We use the github api to fetch latest (up to 100) releases and make a markdown page for them.
"""

import subprocess
from pathlib import Path

import httpx

SAVE_DIR = Path("docs/cli/_auto_generated")
CLI_ENTRY_POINT = Path("packages/divbase-cli/src/divbase_cli/divbase_cli.py")
CLI_COMMANDS_SRC_DIR = Path("packages/divbase-cli/src/divbase_cli/cli_commands")

RELEASES_DIR = Path("docs/releases")
RELEASES_PAGE = RELEASES_DIR / "index.md"
GITHUB_REPO = "ScilifelabDataCentre/divbase"
DIVBASE_RELEASES_URL = f"https://github.com/{GITHUB_REPO}/releases"


SUB_COMMANDS = {
    "auth": CLI_COMMANDS_SRC_DIR / "auth_cli.py",
    "config": CLI_COMMANDS_SRC_DIR / "user_config_cli.py",
    "dimensions": CLI_COMMANDS_SRC_DIR / "dimensions_cli.py",
    "files": CLI_COMMANDS_SRC_DIR / "file_cli.py",
    "query": CLI_COMMANDS_SRC_DIR / "query_cli.py",
    "task-history": CLI_COMMANDS_SRC_DIR / "task_history_cli.py",
    "version": CLI_COMMANDS_SRC_DIR / "version_cli.py",
}


def _fetch_releases() -> list[dict]:
    """Fetch all divbase releases from GitHub API."""
    url = f"https://api.github.com/repos/{GITHUB_REPO}/releases"
    resp = httpx.get(url, headers={"Accept": "application/vnd.github+json"}, timeout=10)
    resp.raise_for_status()
    return resp.json()


# pagination rules of API mean up to 100, I think that is enough for this mirror.
# aka no need to loop through pages here...
_RELEASES_INTRO = f"""# DivBase Releases

Releases are managed on GitHub — the content below is a replica containing up to the 100 latest pulled from the [GitHub releases page]({DIVBASE_RELEASES_URL}) at build time.

"""

_RELEASES_FALLBACK = f"""# DivBase Releases

Releases are managed on GitHub — see the
[GitHub releases page]({DIVBASE_RELEASES_URL}) for the full list.
"""


def generate_releases_page() -> None:
    """Generate docs/releases/index.md from GitHub releases."""
    RELEASES_DIR.mkdir(parents=False, exist_ok=True)
    try:
        releases = _fetch_releases()
    except httpx.HTTPError as exc:
        print(f"WARNING: Could not fetch GitHub releases ({exc}). Writing fallback page.")
        RELEASES_PAGE.write_text(_RELEASES_FALLBACK)
        return

    lines = [_RELEASES_INTRO]
    for release in releases:
        if release.get("draft"):
            continue
        name = release["name"]
        date = release["published_at"][:10]  # YYYY-MM-DD
        url = release["html_url"]
        body = (release.get("body") or "").strip()
        prerelease_label = " *(pre-release)*" if release.get("prerelease") else ""

        lines.append(f"## [{name}]({url}){prerelease_label}\n")
        lines.append(f"*Released: {date} — [View on GitHub]({url})*\n")
        if body:
            lines.append(body + "\n")
        lines.append("---\n")

    print(f"Docs saved to: {RELEASES_PAGE}")
    RELEASES_PAGE.write_text("\n".join(lines))


def generate_cli_docs() -> None:
    """Build CLI documentation using typer's built-in doc generation."""
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

    SAVE_DIR.mkdir(parents=True, exist_ok=True)
    output_file = SAVE_DIR / "divbase-cli.md"
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


def on_startup(command, dirty):
    """
    Build CLI command and release notes documentation using typer's built-in doc generation.

    (The function has to be called this to be picked up by mkdocs as a hook.)

    NOTE: An alternative hook to use could be: 'on_pre_build',
    but that triggers on every save when working on the docs which is a bit excessive.
    """
    generate_releases_page()
    generate_cli_docs()
