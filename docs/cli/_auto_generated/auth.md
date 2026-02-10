# `divbase-cli auth`

Login/logout of DivBase server. To register, visit https://divbase.scilifelab.se/.

**Usage**:

```console
$ divbase-cli auth [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `--install-completion`: Install completion for the current shell.
* `--show-completion`: Show completion for the current shell, to copy it or customize the installation.
* `--help`: Show this message and exit.

**Commands**:

* `login`: Log in to the DivBase server.
* `logout`: Log out of the DivBase server.
* `whoami`: Return information about the currently...

## `divbase-cli auth login`

Log in to the DivBase server.

TODO - think abit more about already logged in validation and UX.
One thing to consider would be use case of very close to refresh token expiry, that could be bad UX.
(But that is dependent on whether we will allow renewal of refresh tokens...)

**Usage**:

```console
$ divbase-cli auth login [OPTIONS] EMAIL
```

**Arguments**:

* `EMAIL`: [required]

**Options**:

* `--password TEXT`: [required]
* `--divbase-url TEXT`: DivBase server URL to connect to.  [default: http://localhost:8000/api]
* `-f, --force`: Force login again even if already logged in
* `--help`: Show this message and exit.

## `divbase-cli auth logout`

Log out of the DivBase server.

**Usage**:

```console
$ divbase-cli auth logout [OPTIONS]
```

**Options**:

* `--help`: Show this message and exit.

## `divbase-cli auth whoami`

Return information about the currently logged-in user.

**Usage**:

```console
$ divbase-cli auth whoami [OPTIONS]
```

**Options**:

* `--help`: Show this message and exit.
