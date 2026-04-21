/* 
Toggle the visibility of the project table based on the selected radio button. 
Used in the personal access token creation page to show/hide project specific table permissions. 
*/
function toggleProjectTable(radio) {
  const table = document.getElementById("projectTable");
  if (table) {
    table.style.display = radio.value === "specific" ? "" : "none";
  }
}
