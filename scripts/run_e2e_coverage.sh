#!/usr/bin/env bash
# Wrapper script to run the DivBase testing suite with coverage.py monitoring distributed across the testing Docker Compose stack
# containers. The script combines the coverage data from the host machine and the containersinto a single HTML report.
#
# The coverage overlay is defined in docker/divbase_compose.tests.coverage.yaml, which extends the base test stack
# to run the FastAPI and Celery worker containers with coverage.py, writing .coverage.* files to a shared docker/coverage-data/ directory.
# The script runs pytest with the --coverage-docker flag to trigger this mode. Since pytest invoke tests that use the CLI to speak to the
# other services in the stack, there will also be a host .coverage file generated for the part of the code that is run on the host machine duriung the test session. 
# After pytest completes, the script copies the host .coverage file into docker/coverage-data/ and runs 'coverage combine' to merge it with the coverage data from the containers. 
# Finally, the script generates an HTML report with paths remapped to the local source tree using the [paths] config in pyproject.toml.
#
# Usage:
#   ./scripts/run_e2e_coverage.sh
#
# Output:
#   htmlcov/index.html   combined HTML coverage report (host + Docker containers)
#   .coverage            merged coverage data file

set -euo pipefail

# Old failed runs may have left coverage data around, so clean it up first to avoid confusing or error-prone results.
echo " - Cleaning up previous coverage data"
rm -f docker/coverage-data/.coverage.*
rm -f .coverage

echo " - Running testing stack (e2e, integration, unit tests) with Docker coverage instrumentation"
PYTEST_EXIT=0
pytest -s tests/ \
    --coverage-docker \
    --cov \
    --cov-branch \
    --cov-context=test \
    --cov-report=term-missing \
    "$@" || PYTEST_EXIT=$?

# Handle how coverage reports should be built given pytest exist codes. Exit codes 0 (all passed) and 1 (some failed) both produce coverage data.
# Codes >=2 are hard errors (interrupted, internal error, no tests collected).
if [[ $PYTEST_EXIT -ge 2 ]]; then
    echo "pytest exited with code $PYTEST_EXIT — skipping coverage report."
    exit $PYTEST_EXIT
fi

# Debug output: if no container coverage data files are found, this will help explain why.
echo ""
echo "- Coverage data files before combine:"
ls -la docker/coverage-data/ || true

echo ""
echo "- Combining host and container coverage data"
cp .coverage docker/coverage-data/.coverage.host
coverage combine docker/coverage-data/

echo " - Generating HTML report"
coverage html --show-contexts

echo ""
echo "Done! Open htmlcov/index.html to browse coverage."
echo "  open htmlcov/index.html"

# Exit with pytest's exit code so test failures are not swallowed by a successful coverage html run.
exit $PYTEST_EXIT
