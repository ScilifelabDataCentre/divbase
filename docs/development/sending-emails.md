# Sending emails

Emails sending for email verification and password resets is handled by the `email_sender.py` service.

Inspiration for the approach was taken heavily from the <https://github.com/fastapi/full-stack-fastapi-template>

## email_verified

- The UserDB Model contains an `email_verified` field and this is checked on login and on refresh token usage.

- The FIRST_ADMIN_USER, test users created by local_dev_setup.py (and equivalent for testing) automatically have their email verified for convenience.

- The admin only api route to create a user can set email_verified to true on creation for convenience.

- If you create a new user by registering via the frontend, this user wont be able to login until they have verified their email.

## email_sender.py

- The service email_sender.py handles sending the different types of emails with emails written using the language: <https://mjml.io/> These are manually converted to Jinja2 templates for actual use in the code.

- Email verification uses JWTs which are included as a link in the sent email. Clicking the link in the received email will call an endpoint with the JWT included as a query param. If token fine, userDB model updated and user can now login.

## MailPit

The docker compose file contains a MailPit service which can be used in local dev/testing (wont be deployed on cluster) for catching all emails sent by DivBase. You can visit the webUI at [http://localhost:8025](http://localhost:8025). Any email sent by DivBase to any email address will be received there so you can use it to check how the emails look and that the verification tokens etc.. work.

## Email verification flow

- Upon account creation, user receives email with link that expires in (currently set to) 24 hours.
- If user visits link within 24 hours, account verified, can log in.
- If user visits after 24 hours, told their link has expired and that they can provide their email to request a new one.
- If a user tries to login who has not verified their email (and they provide correct email + password, they are told they need to verify their email).

## Background tasks

- Email sending is handled as a [background task](https://fastapi.tiangolo.com/tutorial/background-tasks/). This means the response can be returned to the user without waiting for it to be sent.
