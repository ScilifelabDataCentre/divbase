/*  
For both the organisation and role dropdowns in the registration and edit profile forms, 
we want to allow users to select "Other" and then specify their organisation/role if it's not listed.

This handles toggling the new text input for the organisation when "Other" is selected in the organisation dropdown,
and the new text input for the role when "Other" is selected in the role dropdown.

Alongside validation of the inputs on the client side. (Pydantic used for server side validation.)
*/

document.addEventListener("DOMContentLoaded", function () {
  const organisationDropdown = document.getElementById("organisation_dropdown");
  const otherOrganisationDiv = document.getElementById(
    "other_organisation_div",
  );
  const otherOrganisationInput = document.getElementById("organisation_other");

  const roleDropdown = document.getElementById("role_dropdown");
  const otherRoleDiv = document.getElementById("role_other_div");
  const otherRoleInput = document.getElementById("role_other");

  const form =
    document.getElementById("registerForm") ||
    document.getElementById("editProfileForm");
  if (!form) {
    return;
  }

  // Generic validator for "Other" fields
  function validateOtherInput(dropdown, input, message) {
    const value = input.value.trim();
    if (dropdown.value === "Other" && value.length < 3) {
      input.setCustomValidity(message);
    } else {
      input.setCustomValidity("");
    }
  }

  // Specific validator for organisation
  function validateOtherOrganisationInput() {
    validateOtherInput(
      organisationDropdown,
      otherOrganisationInput,
      "Please specify your organisation (at least 3 characters).",
    );
  }

  // Specific validator for role
  function validateOtherRoleInput() {
    validateOtherInput(
      roleDropdown,
      otherRoleInput,
      "Please specify your role (at least 3 characters).",
    );
  }

  // listeners for organisation dropdown
  organisationDropdown.addEventListener("change", function () {
    if (this.value === "Other") {
      otherOrganisationDiv.style.display = "block";
      otherOrganisationInput.required = true;
    } else {
      otherOrganisationDiv.style.display = "none";
      otherOrganisationInput.required = false;
      otherOrganisationInput.value = "";
      otherOrganisationInput.setCustomValidity("");
    }
    validateOtherOrganisationInput();
  });

  otherOrganisationInput.addEventListener("input", function () {
    validateOtherOrganisationInput();
  });

  // listeners for role dropdown
  roleDropdown.addEventListener("change", function () {
    if (this.value === "Other") {
      otherRoleDiv.style.display = "block";
      otherRoleInput.required = true;
    } else {
      otherRoleDiv.style.display = "none";
      otherRoleInput.required = false;
      otherRoleInput.value = "";
      otherRoleInput.setCustomValidity("");
    }
    validateOtherRoleInput();
  });

  otherRoleInput.addEventListener("input", function () {
    validateOtherRoleInput();
  });

  // Form submission, prevent submission if any validation failures.
  form.addEventListener("submit", function (event) {
    validateOtherOrganisationInput();
    validateOtherRoleInput();

    if (!form.checkValidity()) {
      form.reportValidity();
      event.preventDefault();
    }
  });

  // Run on page load to ensure form is in correct state if form submission fails
  // and page/form is re-rendered with previous values.
  // e.g. password don't match etc...
  organisationDropdown.dispatchEvent(new Event("change"));
  roleDropdown.dispatchEvent(new Event("change"));
});
