// styling for options
function addContrastOptions(item, escape) {
  
  
  var clustEl = "<div>" +
    "<div class='input-swatch' style='background-color:" + item.color + "'></div>" +
    escape(item.name) +
    "</div>";
  
  return clustEl;
}



//styling for current item
function addContrastItem(item, escape) {
  
  var clustEl = "<div class='bulk-item'>" +
    "<div class='input-swatch' style='background-color:" + item.color + "'></div>" +
    escape(item.name) +
    "<div class='contrast'>vs</div>" +
    "</div>";
  
  return clustEl;
}



// styling for options
function delContrastOptions(item, escape) {

  // styling if looking at contrast
  var conEl  = "<div title=''>" +
                 "(<div class='input-swatch' style='margin-left: 5px; background-color:" + item.testColor + "'></div>" +
                 " - " +
                 "<div class='input-swatch' style='background-color:" + item.ctrlColor + "'></div>) " +
                 escape(item.test) + " vs " +
                 escape(item.ctrl) +
               "</div>";

  // either cluster or contrast element
  return conEl;
}



//styling for current item
function delContrastItem(item, escape) {

   // styling if looking at contrast
  var conEl  = "<div title=''>" +
                 "(<div class='input-swatch' style='margin-left: 5px; background-color:" + item.testColor + "'></div>" +
                 " - " +
                 "<div class='input-swatch' style='background-color:" + item.ctrlColor + "'></div>) " +
                 escape(item.test) + " vs " +
                 escape(item.ctrl) +
               "</div>";

  // either cluster or contrast element
  return conEl;
}