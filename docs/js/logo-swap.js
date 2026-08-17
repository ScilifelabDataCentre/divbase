/*
Swaps the logo from dark to light mode based on user preference.
see related: docs/overrides/partials/logo.html
*/
const LIGHT_LOGO = "assets/img/scilifelab_logo_light.svg";
const DARK_LOGO = "assets/img/scilifelab_logo_dark.webp";

function updateLogo() {
  const logo = document.querySelector('a[data-md-component="logo"] img');
  if (!logo) {
    return;
  }

  const targetFile =
    document.body?.getAttribute("data-md-color-scheme") === "slate"
      ? DARK_LOGO
      : LIGHT_LOGO;

  logo.src = new URL(targetFile, __md_scope).toString();
}

const observer = new MutationObserver(updateLogo);

observer.observe(document.body, {
  attributes: true,
  attributeFilter: ["data-md-color-scheme"],
});

updateLogo();

if (typeof document$ !== "undefined") {
  document$.subscribe(updateLogo);
} else {
  document.addEventListener("DOMContentLoaded", updateLogo, { once: true });
}
