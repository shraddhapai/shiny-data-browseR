$(document).on('click', '.selectable div table tbody tr', function(e){
	var el = $(this);
	if (!e.shiftKey){
		$(this).siblings().removeClass("rowsSelected");
	}
	$(this).addClass("rowsSelected", this.clicked);
	el.trigger("change");
});	

//$(document).on('load', 'table', function(e) {
//var table=document.getElementById("o_sampleTable");
//for (var i=0,row; row=table.rows[i]; i++) {
//	$row.addClass("rowsSelected");
//}
//table.trigger("change");
//}

var selectRowBinding = new Shiny.InputBinding();
$.extend(selectRowBinding, { // SP: methods associated with selected rows.
	find: function(scope) {
		return $(scope).find(".selectable");
	},
	getValue: function(el){	// SP: print table method.
    tbl = $(el).find("table");
    var out = [];
    $rows = $(tbl).children().children('.rowsSelected');
    if($rows.length == 0) return -1; //SP: Return -1 if no rows selected

//    uncomment to return row numbers instead 
//    (would be wrong if table is sorted or filtered)

//    $rows.each(function(row, v){
//      out[row] = $(v).index() + 1;
//    });
//  return out;    

	// SP:iterate over rows
    $rows.each(function(row,v) {
      $(this).find("td").each(function(cell,v) {
        if (typeof out[row] === 'undefined') out[row] = [];
        out[row][cell] = $(this).text(); 
      });
    });
    return out; //SP: responsible for printing selected rows.
	},
	setValue: function(el, value) {
	},
	subscribe: function(el, callback) {
		$(el).on("change.selectRowBinding", function(e) {
			callback();
		});
	},
	unsubscribe: function(el) {
	  $(el).off(".selectRowBinding");
	}
});
Shiny.inputBindings.register(selectRowBinding);
