/*
This script provides a way to toggle the theme (light/dark mode) once the page is loaded.
This script works with "determine_theme.js" which runs immeditaly and sets the initial user prefered theme in local storage (if not already there).

The html associated with this JS is in a Hugo partial "change_theme.html", which is used in the navbar.
*/
document.addEventListener("DOMContentLoaded", () => {
  const htmlElement = document.documentElement;
  const themeOptions = document.querySelectorAll(
    ".dropdown-item[data-bs-theme-value]",
  );

  // helper function to swap image srcs.
  const updateImages = (theme) => {
    let srcAttr = `data-${theme}-src`;
    document.querySelectorAll("img[data-dark-src]").forEach((img) => {
      img.src = img.getAttribute(srcAttr);
    });
  };

  // helper function to update theme dropdown with which color theme is currently active
  const setActiveClass = (theme) => {
    themeOptions.forEach((option) => {
      if (option.getAttribute("data-bs-theme-value") === theme) {
        option.classList.add("active");
      } else {
        option.classList.remove("active");
      }
    });
  };

  const updateTheme = (theme) => {
    htmlElement.setAttribute("data-bs-theme", theme);
    localStorage.setItem("bsTheme", theme);
    updateImages(theme);
    setActiveClass(theme);
  };

  // Initial setup
  const currentTheme = htmlElement.getAttribute("data-bs-theme");
  setActiveClass(currentTheme);
  updateImages(currentTheme);

  // Event listener for theme change
  themeOptions.forEach((item) => {
    item.addEventListener("click", () => {
      const newTheme = item.getAttribute("data-bs-theme-value");
      updateTheme(newTheme);
    });
  });
});
