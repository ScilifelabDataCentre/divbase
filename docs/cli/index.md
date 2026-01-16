# DivBase CLI Command Reference

This section of the documentation covers all available CLI commands for `divbase-cli`.

Each page in this section of the documentation is auto-generated from the divbase-cli code and is also available directly via the CLI tool by running:

```bash
divbase-cli --help # or -h
```

You can get more specific help for a command by specifying any divbase command and `--help` or `-h`. For example:

```bash
divbase-cli config -h # help on all divbase-cli config commands
divbase-cli config add -h # help for the divbase-cli config add command.
```

`divbase-cli` is made up of a series of subcommands grouped into categories. These categories are:

- `auth`: Commands for logging in and out of DivBase.
- `config`: Manage your user configuration file. Your configuration file contains information about which DivBase projects you're a member of.
- `dimensions`: Create and inspect dimensions (number of samples, number of variants, scaffold names) of the VCF files in a project"
- `files`: Download/upload/list files to/from the project's store on DivBase.
- `queries`: Run queries on the VCF files stored in the project's data store on DivBase. Queries are run on the DivBase server.
- `task-history`: Get the task history of query jobs submitted by the user to the DivBase server.
- `versions`: Add, view and remove versions representing the state of all files in the entire project at the current timestamp.
