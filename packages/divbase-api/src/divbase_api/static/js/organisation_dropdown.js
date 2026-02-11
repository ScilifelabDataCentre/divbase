/*  
Toggles display of a text input for a user to specify their organisation if 
they choose the option "Other" from the dropdown of swedish universities.
Used in both registration and edit profile pages.
*/

document.addEventListener("DOMContentLoaded", function () {
  const organisationDropdown = document.getElementById("organisation_dropdown");
  const otherOrganisationDiv = document.getElementById("other_organisation_div");
  const otherOrganisationInput = document.getElementById("organisation_other");

  // Check if the elements exist before adding listeners
  if (organisationDropdown && otherOrganisationDiv && otherOrganisationInput) {
    organisationDropdown.addEventListener("change", function () {
      if (this.value === "Other") {
        otherOrganisationDiv.style.display = "block";
        otherOrganisationInput.required = true;
      } else {
        otherOrganisationDiv.style.display = "none";
        otherOrganisationInput.required = false;
        otherOrganisationInput.value = ""; // Clears the value if hidden
      }
    });
    // Run on page load to ensure form is in correct state if form submission fails 
    // and page/form is re-rendered with previous values.
    // e.g. password don't match etc... 
    organisationDropdown.dispatchEvent(new Event("change"));
  }
});