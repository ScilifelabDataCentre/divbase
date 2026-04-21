# Tips for using DivBase programmatically

Below is a set of tips for users who want to use DivBase in scripts/pipelines/programmatically.

If you have any tips or suggestions to add to this page or any desired features, please let us know!

TODO: E.G. how to wait for a query to be complete, and download the query results programmatically.

## Parse divbase-cli files ls/info output programmatically

1. You can make the output of the `divbase-cli files info` and `divbase-cli files ls` commands in TSV format for easier parsing. Use the `--tsv` flag:

    ```bash
    divbase-cli files ls --tsv
    divbase-cli files info FILE_NAME --tsv
    ```

    You can do the same for any [project versions](./project-versioning.md) you've created for your project:

    ```bash
    divbase-cli version ls --tsv
    divbase-cli version info VERSION_NAME --tsv
    ```

2. Rather than first downloading a file, you can stream a file from the command line and pipe it into other tools for processing directly without saving it to disk.

    ```bash
    divbase-cli files stream my_file.vcf.gz | zcat | less
    ```

    !!! Info
        BCFTools accepts stdin as input, so you can also pipe a VCF file directly into BCFTools without saving it first:

        ```bash
        divbase-cli files stream my_file.vcf.gz | bcftools view -h -
        ```

## Login programmatically

If you want to use DivBase in a script/job that may take **longer than 1 week** to complete from the time of submission, you can login to DivBase programmatically in the script itself. For security reasons (not having to store your password) it is preferable to instead log in interactively before submitting the script, which will keep you logged in for 1 week.

To login programmatically, use the `--password-stdin` (or `-p`) flag to provide your password via standard input (STDIN) when logging into DivBase with the `divbase-cli auth login` command. Using STDIN prevents the password from ending up in the shell's history, or log-files.

If you use a password manager that allows you to output your password to standard output, you can pipe the password directly into the `divbase-cli auth login` command. For example, if you use `pass` as your password manager, you could do:

```bash
pass show my_divbase_password | divbase-cli auth login EMAIL_ADDRESS --password-stdin
```

An alternative (but less secure) way to do this is to store your password in a plain text file (make sure to set the appropriate permissions on the file to help make it more secure) and then read the password from the file and pipe it into the `divbase-cli auth login` command. For example:

```bash
chmod 600 ~/my_password.txt  # only readable/writeable by the owner
cat ~/my_password.txt | divbase-cli auth login EMAIL_ADDRESS --password-stdin
```

Another less secure alternative is to use environment variables to store your password and then echo the password from the environment variable and pipe it into the `divbase-cli auth login` command. For example:

```bash
echo $DIVBASE_PASSWORD | divbase-cli auth login EMAIL_ADDRESS --password-stdin
```

Where `DIVBASE_PASSWORD` is an environment variable that contains your password.
