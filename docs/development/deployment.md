# Deployment of DivBase Server

DivBase is deployed on a kubernetes cluster and its deployment is managed in our [private repository, argocd-divbase](https://github.com/ScilifelabDataCentre/argocd-divbase), see there for details.

In this repository there are 3 things related to deployment:

1. Release the docker images of the api and worker to the GitHub Container Registry (GHCR)
2. Publish new versions of the divbase-cli and divbase-lib packages to PyPI.
3. Publish the documentation to GitHub Pages (which is used for hosting the docs site).

For releases, we have a GitHub Actions workflow (`release.yaml`) that does all three of these things automatically when the release is published.

- For development and testing purposes, you can run specific github actions via `workflow_dispatch`:

  - `publish_dev_images.yaml`: Build and push new docker images to the GitHub Container Registry. Images will be tagged with commit hash and latest-dev. A pytest run will run after the images are built to test the new images.
  - `publish-to-test-pypi.yaml`: Publish divbase-lib and divbase-cli packages to **Test** PyPI.
  - `publish_docs.yaml`: Publish a new version of the docs site (we have no dev docs site, just 1 single shared docs site).

- To trigger a workflow_dispatch, head to the "Actions tab" on GitHub and to the action name of choice.

## Publishing of `divbase-cli` and `divbase-lib`

The packages `divbase-cli` and `divbase-lib` are published to PyPI for user installation via pip/uv etc...

We have accounts on both [PyPI](https://pypi.org/) and [TestPyPI](https://test.pypi.org). With TestPyPI you can publish test versions of the packages before pushing to the main PyPI repository. Users wont accidentally install versions on the test instance unless they explicitly specify to install from there.

- [divbase-cli on PyPI](https://pypi.org/project/divbase-cli/)
- [divbase-lib on PyPI](https://pypi.org/project/divbase-lib/)
- [divbase-cli on test PyPI](https://test.pypi.org/project/divbase-cli/)
- [divbase-lib on test PyPI](https://test.pypi.org/project/divbase-lib/)

(Login credentials are in bitwarden, but if you use the GH actions below you don't need to use them.)

### Good to know

1. You cannot upload the same version of a package to PyPI or TestPyPI if that version already exists (even if you delete the existing package). So in such cases you need to bump the version number in the packages before uploading new versions.

2. We set the version numbers of `divbase-{lib/cli/api}` to be identical for any release. This is done for the sake of simplicity (obvious which version of which package supports each other). So when you bump a version number for a new release, make sure to bump it in all packages to the same value, (read below for a helper script that does this for you).

### Publishing to TestPyPI

To test how the installation of a new version would be without affecting the production PyPI registry:

1. Bump the version number in each packages `__init__.py` file and the version of `divbase-lib` to install in the `pyproject.toml` for divbase-cli. You can do this easily with the helper script:

    ```bash
    uv run scripts/set_divbase_version.py [NEW_VERSION]
    ```

    !!! Tip "How should I bump the version when I'm just testing something?"
        If a semantic version bump does not make sense (perhaps you're just testing something not ready for end users), then you can append a dev suffix to the current version, e.g. if the current version is `1.2.3`, you can set the new version to `1.2.3.dev1` or `1.2.3.dev2` etc...

2. Validate you bumped versions correctly across all files by running:

    ```bash
    uv run scripts/ci/validate_package_versions.py
    ```

    (gh actions will also validate this when you trigger the workflow.)

3. Commit and push your changes to GitHub (can be to main or a feature branch).

4. Go the actions tab and run the workflow **Publish divbase-lib and divbase-cli packages to Test PyPI**.

    !!! Warning
        Make sure to run the one to the TEST PyPI, not to production PyPI.

5. Once the action has run, you can now test the installation of the new version.

    !!! Warning "Make sure you're testing the installed version"
        To make sure you're testing the installed version from TestPyPI, and not the version from your local env in your local git repo, you should create a new terminal window (not in your divbase git repo root) and run the install command there.

    ```bash
    # open a new terminal and not the in the root of the divbase git repo
    uv tool install --index https://test.pypi.org/simple --index-strategy unsafe-best-match divbase-cli==[version]

    # run your divbase commands now...
    divbase --version
    ```

    When finished testing, you can uninstall the test version with:

    ```bash
    uv tool uninstall divbase-cli
    ```

### Publishing to PyPI

The production release is automated, once a new GH release is published on GitHub, the actions workflow (`release.yaml`) will run to publish the packages to PyPI.

If something goes wrong with the GH action that automatically publishes to PyPI:

1. Did the upload to PyPI succeed? If so, we can consider "yanking" the release on PyPI to prevent users from installing a broken version while we fix the issue.
2. Commit the fix, and bump the version number accordingly.
3. Make a new release with the fix, this will trigger a new GH action to publish the new version to PyPI.

!!! Info
    If really needed, you can consider enabling `workflow_dispatch` on the GH action workflow to manually trigger it.
    Or locally run:
    ```bash
    # Make sure to have the correct version numbers set in the packages before running these commands

    uv build --package divbase-cli --no-sources
    uv build --package divbase-lib --no-sources
    uv publish .....
    ```
    [see here for why you need to add the flag --no-sources to the build steps](https://docs.astral.sh/uv/guides/package/#building-your-package).
    But this should really not be needed in normal cases.

## Considerations for Kubernetes deployment

TODO: add some lessons learned from the k8s deployment, e.g. on resources needed for the celery worker.
