---
search:
  exclude: true
---

# Documentation Generation

Documentation for DivBase is generated using [Zensical](https://zensical.org/).

The `zensical.toml` at the root of the repository contains the documentation configuration for the site. All content inside the `docs/` folder can be used to create documentation pages.

To add a new page to the documentation site you must add it to the `nav` section of the `zensical.toml` file.

To make the documentation more user friendly we don't include the developer specific section of the docs in the site search. To do this add:

```
---
search:
  exclude: true
---
```

To any new developer documentation file.

## View documentation locally

To serve the documentation locally (with hot reloading):

```bash
uv run scripts/build_docs.py && uv run zensical serve
```

The docs site will be available at: <http://localhost:8008/divbase/>

If you make changes to the documentation files inside the `docs/` folder, the local server will automatically reload to reflect those changes. Changes to the `zensical.toml` file require a restart of the server to take effect.

!!! info "You technically don't need to run the build script first in local development"

    - The `scripts/build_docs.py` builds some auto-generated pages (`divbase-cli` command reference and a page of recent GitHub releases).
    - These pages will just return a 404 if you're try to access them locally if you haven't run the build script.
    - The auto-generated pages are .gitignore'd, and the GH action run to deploy the docs always builds them before a new release (so they are always in sync).

## Deploying the documentation site

A GitHub Actions workflow is set up to build and deploy the documentation site to GitHub Pages on:

1. A new release `.github/workflows/release.yaml`.
2. A workflow dispatch (trigger manually in gh actions tab for any branch you like). `.github/workflows/publish_docs.yaml`.
