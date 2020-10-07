function bulkContrastOptions(item, escape) {
  
  
  var clustEl = "<div>" +
    "<div class='input-swatch' style='background-color:" + item.color + "'></div>" +
    escape(item.name) +
    "</div>";
  
  return clustEl;
}



//styling for current item
function bulkContrastItem(item, escape) {
  
  var clustEl = "<div class='bulk-item'>" +
    "<div class='input-swatch' style='background-color:" + item.color + "'></div>" +
    escape(item.name) +
    "<div class='contrast'>vs</div>" +
    "</div>";
  
  return clustEl;
}
