# Account Management

## Getting a project

!!! Warning
    DivBase is under active development and not in production yet. We are shortly planning to create a number of projects for some pilot users. If you would like to be involved in this pilot testing phase, please let us know!

Projects in DivBase are created by the DivBase team. To request a new project, you need to contact us. Please include a very short description of what the project would be for and who you plan to give access too. Please email us with your university affliated email address. please make sure you have already registered for a DivBase account and verified your email before contacting us about creating a project.

Once your project is created, you will be added as a member with the **Manage** role, which allows you to add and remove other members.

## Adding members to a project

Project members can be added by anyone with the **Manage** role on that project.

1. Navigate to your project page via **Projects** in the navigation bar.
2. Open the project you want to manage.
3. In the **Members** section, click **Add Member**.
4. Enter the email address of the person you want to add and select their role (`read`, `edit`, or `manage`).
5. Click **Add Member** to confirm.

!!! note "The user must already have a DivBase account"
    The person must have registered and verified their email before they can be added. If they haven't signed up yet, ask them to do so first.

You can also change a member's role or remove them from the project from the same Members section.

## Resetting your password

If you've forgotten your password, click **Forgot your password?** on the login page. Enter your email address and you'll receive a reset link. The link expires after a short period — if it has expired you can simply request a new one.

You can also change your password at any time from your **Profile** page once logged in.

## Personal Access Tokens

Personal Access Tokens (PATs) let you authenticate with DivBase programmatically without having to store your password on the device. This could be useful in pipelines, or HPC jobs that need to interact with DivBase. You can manage your personal access tokens from the **Profile** → **Personal Access Tokens** page.

!!! Info
    To learn more about how to use a PAT in scripts/pipelines and/or HPC jobs, see [Using DivBase Programmatically](./using-divbase-programmatically.md#use-personal-access-tokens-to-authenticate-programmatically).

To create a PAT:

1. Log in to DivBase and navigate to your **Profile** page. From their you'll see a section to manage your Personal Access Tokens.
2. When adding a new token, you need to give it a name, an optional description and an expiry date.
3. You can also configure the scope of the PAT (What permissions the token has):
    - **All projects** — the token can access all projects you're a member of, with same role as you have.
    - **Specific projects** — restrict the token to selected projects, each with a maximum role.
    - **Task history** — Allow you to use the token to view your entire task history, not filtered by project.
4. Click **Create**. Copy the token immediately — it is shown **only once** and cannot be recovered afterwards.

You can have up to **5 active tokens** at a time. To revoke a token, click **Revoke** next to it on the Personal Access Tokens page.

For how to use a PAT in scripts and HPC jobs, see [Using DivBase Programmatically](./using-divbase-programmatically.md#use-personal-access-tokens-to-authenticate-programmatically).
