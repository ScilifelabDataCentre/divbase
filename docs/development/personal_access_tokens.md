# Personal Access Tokens (PATs)

Personal Access Tokens (PATs) provide an alternative authentication mechanism to password and JSON web tokens (JWTs) for users who need to make use of DivBase programmatically. For example, from HPC job scripts or pipelines.

A user never logs in with a PAT, they instead create a PAT via the frontend and pass the PAT in each API call. PATs can be configured to expire or not expire. On every API call, a db lookup is made in the `personal_access_tokens` table to verify the token and check its permissions - and get the corresponding user. This means that PATs can be revoked immediately by deleting/soft-deleting them from the db.

## PATs vs JWTs

JWTs: When an interactive user logs in via the CLI or the frontend, the server issues a short-lived **JWT access token** and a longer-lived **JWT refresh token**. The CLI auto-refreshes the access token when needed using the refresh token.

PATs have no refresh cycle and no concept of login. A PAT is much closed to an API key. A user creates a PAT via the frontend, copies the plaintext token to their device (probably making it an environment variable) and each API call to DivBase will use the PAT as a `Bearer` token in the `Authorization` header. The `get_current_user` in `deps.py` is responsible for handling the two diff auth mechanisms in the API (JWTs or PATs). All PATs are prefixed with `divbase_pat_` followed by a long random string. The token is only shown to the user at creation time and a hash of it is stored in the db.

## Permissions model

Permissions are stored in the `permissions` column of the `personal_access_tokens` table as a **JSONB**, validated against the `PATPermissions` Pydantic model (defined in `divbase-lib`):

```python
class PATPermissions(BaseModel):
    all_projects: bool = False
    projects: dict[str, str] = Field(default_factory=dict)
    task_history: bool = False
    whoami: bool = False
```

The JSONB is nullable, and a `NULL` column means same access level as the user.

| Scope | Meaning |
|---|---|
| `all_projects` | Access all projects the user is a member of, at their membership role |
| `projects` | Access specific projects only; each project mapped to a max role (`read`, `edit`, `manage`) |
| `task_history` | Access the task history endpoints which are not project specific |
| `whoami` | Always `True` — granted implicitly on all PATs, not user-configurable |

The effective project role for a scoped PAT is `min(pat_role, user_membership_role)` — the PAT can never escalate beyond the user's actual membership.
