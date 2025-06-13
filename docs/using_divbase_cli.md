
# Using divBase CLI with the dev cluster deployment of MinIO

Make sure you have followed the install instructions in the main README.md to install the divbase CLI before continuing.

## 1. Login to your account and create and download your access and secret keys

With the credentials you recieved from someone in the team, login to the MinIO instance at:

<https://minio.divbase-testground.scilifelab-2-dev.sys.kth.se/browser>

Go to Access Keys (left side menu) and click create a new access key and secret key. Follow the default settings, and make sure to copy the access key and secret key.

These need to be available as environment variables in order to use the CLI tool.

(One way to do this) In the root of the repository, create a file called `.env` and add the following lines:

```bash
DIVBASE_ACCESS_KEY=[YOUR ACCESS KEY ADDED HERE, WITHOUT BRACKETS]
DIVBASE_SECRET_KEY=[YOUR SECRET KEY ADDED HERE, WITHOUT BRACKETS]
```

This `.env` file is personal and *should not* be committed to the repository (it is `.gitignored` already). This `.env` file will be automatically used by the CLI tool to give you access to the bucket.

## 2. First time running the CLI Tool

Each CLI command comes with a `--help` or `-h` option.

To see the available commands and options, run:

```bash
python -m divbase-cli -h
```

### 2.1. Start by creating a user configuration file. This is done by running

```bash
python -m divbase-cli config create
```

Which will create a file in your home dir with path: `~/.divbase-cli/config.yaml`

To see the contents of this file, either go look at it or run:

```bash
python -m divbase-cli config show
```

### 2.2. Add a default bucket to the configuration file

This isn't essential, but means you don't have to specify the bucket name every time you run a command.

```bash
python -m divbase-cli config add-bucket [BUCKET_NAME] --default
```

By adding the `--default` or `-d` flag, this bucket will be set as the default bucket for all commands (you can have multiple buckets)

If you run:

```bash
python -m divbase-cli config show
```

You'll now see your bucket added to the configuration file and marked as the default bucket.

### 2.3. If you're in a new project, you can create the bucket versioning file by running

```bash
python -m divbase-cli version create
```

This will create a file in the bucket called `.bucket_versions.yaml`

An already existing bucket likely has this file, we can see the contents of this file by running:

```bash
python -m divbase-cli version list
```

These are user specified version of the entire buckets state at a given time. This is useful for keeping track of changes to the bucket or marking important states of the bucket. You do not need to specify this if you're only interested in working with the latest versions of the files in the bucket.

If after working with the bucket for a while you want to version the current state of the bucket, you can run:

```bash
python -m divbase-cli version add [OPTIONS] NAME
```

### 2.4. To upload and download files to and from the bucket, you can use the `file` subcommand

```bash
python -m divbase-cli file -h
```

Remember that unless you specified the bucket you want to use in the command, the default bucket set in your user config will be used.

*Bonus:*

To download files from bucket at a specific bucket version/state, we can use the --bucket-version option and specify the version name we want to download from:

```bash
python -m divbase-cli file download file1.vcf.gz file2.vcf.gz --bucket-version=v0.1.0
```
