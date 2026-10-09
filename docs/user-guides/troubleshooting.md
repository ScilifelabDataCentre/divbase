# Troubleshooting

This page collects issues that users might encounter when using DivBase.

Please feel free to contact us at [datacentre@scilifelab.se](mailto:datacentre@scilifelab.se) if you're struggling!

## User accounts and login

??? question "Authentication Issues"
    - Make sure you've verified your email address
    - Check if you can login to your account on the [DivBase Website](https://divbase.scilifelab-2-prod.sys.kth.se). If it fails on the website it will also fail on the CLI.
    - Try logging out and back in: `divbase-cli auth logout` then `divbase-cli auth login`

## Jobs

??? question "My job exceeded the maximum run time"
    Jobs are stopped when they reach their [maximum run time](running-queries-overview.md#maximum-run-times). Try splitting a large query into smaller ones, or contact us if you need to run longer jobs.

??? question "The worker running this job stopped unexpectedly"
    This means an unexpected error caused the job to stop mid-run. This is not expected to happen often and does not imply something wrong with your job.
    Please resubmit it. If it happens repeatedly, please contact us.

## Debugging `divbase-cli` errors

By default expected `divbase-cli` errors are printed in a user-friendly format without including the full traceback.

If you want to see the full traceback for debugging purposes, you can set the environment variable `DIVBASE_TRACEBACKS_ON=1` before running your command.

For example, this command will fail and show the full traceback:

```bash
DIVBASE_TRACEBACKS_ON=1 divbase-cli files mkdir 'dir-?with:invalid:ch@rs'
```
