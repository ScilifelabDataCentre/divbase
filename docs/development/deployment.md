# Deployment of DivBase Server

- DivBase is deployed on a kubernetes cluster and its deployment is managed in our [private repository, argocd-divbase](https://github.com/ScilifelabDataCentre/argocd-divbase), see there for details.

- In this repository the main thing that needs to be done is to build and push the new docker images for the API and worker.

- We have a GitHub Actions workflow that does this (`.github/workflows/publish_images.yaml`). You can trigger a workflow_dispatch to manually create the docker images at a given commit.

- Head to the actions tag on GitHub and to the action "Build and push both docker images to the GitHub Container Registry". From there click run manual workflow. You can choose to specify the name of the image tag if you want. Otherwise leave the input blank and it will be tagged with the full commit hash.

## Deployment of DivBase CLI

The packages divbase-cli and divbase-lib are published to PyPI for user installation via pip/uv etc...

We have accounts on both [PyPI](https://pypi.org/) and [TestPyPI](https://test.pypi.org). With TestPyPI you can publish test versions of the packages before pushing to the main PyPI repository. Users wont accidentally install versions on the test instance unless they explicitly specify to install from there.

- [divbase-cli on PyPI](https://pypi.org/project/divbase-cli/)
- [divbase-lib on PyPI](https://pypi.org/project/divbase-lib/)
- [divbase-cli on test PyPI](https://test.pypi.org/project/divbase-cli/)
- [divbase-lib on test PyPI](https://test.pypi.org/project/divbase-lib/)

(Login credentials are in bitwarden, but if you use the GH actions below you don't need to use them.)

Something good to know is that you cannot upload the same version of a package to PyPI or TestPyPI if that version already exists (even if you delete the existing package). So in such cases you need to bump the version number in the packages before uploading new versions.

### Publishing to TestPyPI

To test how the installation of a new version would be without affecting the production PyPI registry:

1. Bump the version number in each packages `__init__.py` file and the version of `divbase-lib` to install in the `pyproject.toml` for divbase-cli. So these 4 files:
    - `packages/divbase-api/src/divbase_api/__init__.py`
    - `packages/divbase-lib/src/divbase_lib/__init__.py`
    - `packages/divbase-cli/src/divbase_cli/__init__.py`
    - `packages/divbase-cli/pyproject.toml` (the version of `divbase-lib` to install)

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
        To make sure you're testing the installed version from TestPyPI, and not the version in your local git repo, you should create a new terminal window (not in your divbase git repo virtual env) and run the install command there.

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

The production release is automated, once a new release is published on GitHub, a [GitHub Actions workflow will run to publish the packages to PyPI](.github/workflows/publish-to-pypi.yaml).

If something goes wrong with the GH action that automatically publishes to PyPI:

1. Did the upload to PyPI succeed? If so, we can consider "yanking" the release on PyPI to prevent users from installing a broken version while we fix the issue.
2. Commit the fix, and bump the version number accordingly.
3. Make a new release with the fix, this will trigger a new GH action to publish the new version to PyPI.

!!! Info
    If really needed, you can consider enabling `workflow_dispatch` on the GH action workflow to manually trigger it.
    Or locally use `uv build` and then `uv publish` to manually build and publish the packages to PyPI.
    But this should really not be needed in normal cases.
