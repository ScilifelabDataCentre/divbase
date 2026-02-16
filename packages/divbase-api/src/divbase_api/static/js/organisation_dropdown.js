/*  
Toggles display of a text input for a user to specify their organisation if 
they choose the option "Other" from the dropdown of swedish universities.
Used in both registration and edit profile pages.
*/

document.addEventListener("DOMContentLoaded", function () {
  const organisationDropdown = document.getElementById("organisation_dropdown");
  const otherOrganisationDiv = document.getElementById("other_organisation_div");
  const otherOrganisationInput = document.getElementById("organisation_other");

  const form = document.getElementById("registerForm") || document.getElementById("editProfileForm");
  if (!form) {
    return;
  }

  function validateOtherOrganisation() {
    const otherOrganisation = otherOrganisationInput.value.trim();
    if (
      organisationDropdown.value === "Other" &&
      (!otherOrganisation || otherOrganisation.length < 3)
    ) {
      otherOrganisationInput.setCustomValidity("Please specify your organisation (at least 3 characters).");
    } else {
      otherOrganisationInput.setCustomValidity("");
    }
  }

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
  });

  otherOrganisationInput.addEventListener("input", validateOtherOrganisation);

  // Validation for "Other" organisation input on form submission will prevent form submission
  form.addEventListener("submit", function (event) {
    validateOtherOrganisation();

    if (!form.checkValidity()) {
      form.reportValidity();
      event.preventDefault();
    }
  });

  // Run on page load to ensure form is in correct state if form submission fails 
  // and page/form is re-rendered with previous values.
  // e.g. password don't match etc... 
  organisationDropdown.dispatchEvent(new Event("change"));
});