#!/usr/bin/env bash
# Wrapper script to run e2e tests with coverage.py distributed across the testing Docker Compose stack
# containers and to combine the coverage data into a single HTML report.
#
# The coverage overlay is defined in docker/divbase_compose.tests.coverage.yaml, which extends the base test stack
# to run the FastAPI and Celery worker containers with coverage.py, writing .coverage.* files to a shared docker/coverage-data/ directory.
# The script runs pytest with the --coverage-docker flag to trigger this mode. Since pytest invoke tests that use the CLI to speak to the
# other services in the stack, there will also be a host .coverage file generated for the part of the code that is run on the host machine duriung the test session. 
# After pytest completes, the script copies the host .coverage file into docker/coverage-data/ and runs 'coverage combine' to merge it with the coverage data from the containers. 
# Finally, the script generates an HTML report with paths remapped to the local source tree using the [paths] config in pyproject.toml.
#
# Usage:
#   ./scripts/run_e2e_coverage.sh                      # all e2e tests (excludes playwright)
#   ./scripts/run_e2e_coverage.sh -k test_auth_cli     # filter to specific tests
#   ./scripts/run_e2e_coverage.sh --run-slow           # include slow tests
#
# Output:
#   htmlcov/index.html   combined HTML coverage report (host + Docker containers)
#   .coverage            merged coverage data file

set -euo pipefail

echo " - Running e2e tests with Docker coverage instrumentation"
PYTEST_EXIT=0
pytest tests/e2e_integration/ \
    --ignore=tests/e2e_integration/playwright \
    --coverage-docker \
    --cov \
    --cov-context=test \
    --cov-report=term-missing \
    "$@" || PYTEST_EXIT=$?

# Handle how coverage reports should be built given pytest exist codes. Exit codes 0 (all passed) and 1 (some failed) both produce coverage data.
# Codes >=2 are hard errors (interrupted, internal error, no tests collected).
if [[ $PYTEST_EXIT -ge 2 ]]; then
    echo "pytest exited with code $PYTEST_EXIT — skipping coverage report."
    exit $PYTEST_EXIT
fi

echo ""
echo "- Combining host and container coverage data"
cp .coverage docker/coverage-data/.coverage.host
coverage combine docker/coverage-data/

echo " - Generating HTML report"
coverage html

echo ""
echo "Done! Open htmlcov/index.html to browse coverage."
echo "  open htmlcov/index.html"

# Re-exit with pytest's exit code so callers (CI) see test failures.
exit $PYTEST_EXIT
