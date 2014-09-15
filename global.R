#selDataTableDef <- function() {
#cat("Calling this")
#}

selDataTableOutput <- function (outputId) 
{
  x <- tagList(singleton(shiny::tags$head(
  	shiny::tags$link(rel = "stylesheet", 
    	type = "text/css", href = "shared/datatables/css/DT_bootstrap.css"),
    shiny::tags$style(type="text/css", ".rowsSelected td{
               background-color: rgba(131,112,225,0.2) !important}"),
    shiny::tags$style(type="text/css", ".selectable div table tbody tr{
               cursor: hand; cursor: pointer;}"),
    shiny::tags$style(type="text/css",".selectable div table tbody tr td{
               -webkit-touch-callout: none;
               -webkit-user-select: none;
               -khtml-user-select: none;
               -moz-user-select: none;
               -ms-user-select: none;
               user-select: none;}"),                          
	shiny::tags$style(HTML("
html {
  font-family: sans-serif;
  -ms-text-size-adjust: 100%;
  -webkit-text-size-adjust: 100%;
}
						   
body {
  margin: 0px;
  margin-top:80px;
  margin-left:10px;
  width:1280px;	
}

.header-fixed {
	width:100%;
	position:fixed;
	top:0px
}
.header {
	height:70px;
	border:1px solid #CCC;
	width:1280px;
	margin: 0px auto;
	z-index:100;
	opacity:0.8;
}

.pagebar-col {
  background-color: #458cc3;
  border-color: #727272;
  color: #eeeeee;
  font-weight:800;
  border:3px solid transparent;
}
.pagebar-col2 {
  background-color: #2c5B80;
  border-color: #727272;
  color: #eeeeee;
  font-weight:400;
  border:3px solid transparent;
}


.btn {
  display: inline-block;
  margin-bottom: 0;
  font-weight: normal;
  text-align: center;
  vertical-align: middle;
  cursor: pointer;
  background-image: none;
  border: 3px solid transparent;
  white-space: nowrap;
  padding: 6px 12px;
  font-size: 14px;
  line-height: 1.42857143;
  border-radius: 10px;
  -webkit-user-select: none;
  -moz-user-select: none;
  -ms-user-select: none;
  user-select: none;
}

.btn-success {
  color: #ffffff;
  background-color: #28b62c;
  border-color: #23a127;
}

.btn-success:hover,
.btn-success:focus,
.btn-success:active,
.btn-success.active,
.open > .dropdown-toggle.btn-success {
  color: #ffffff;
  background-color: #1f8c22;
  border-color: #186f1b;
}
.btn-success:active,
.btn-success.active,
.open > .dropdown-toggle.btn-success {
  background-image: none;
}
.btn-success.disabled,
.btn-success[disabled],
fieldset[disabled] .btn-success,
.btn-success.disabled:hover,
.btn-success[disabled]:hover,
fieldset[disabled] .btn-success:hover,
.btn-success.disabled:focus,
.btn-success[disabled]:focus,
fieldset[disabled] .btn-success:focus,
.btn-success.disabled:active,
.btn-success[disabled]:active,
fieldset[disabled] .btn-success:active,
.btn-success.disabled.active,
.btn-success[disabled].active,
fieldset[disabled] .btn-success.active {
  background-color: #28b62c;
  border-color: #23a127;
}
.btn-success .badge {
  color: #28b62c;
  background-color: #ffffff;
}

.btn-lg {
  padding: 10px 16px;
  font-size: 16px;
  line-height: 1.33;
  border-radius: 6px;
  font-weight:600;
	
}

.btn-xlg {
  padding: 10px 16px;
  font-size: 20px;
  line-height: 1.33;
  border-radius: 10px;
	
}
.text-primary,
.text-primary:hover {
  color: #158cba;
}

.panel-footer {
  padding: 10px 15px;
  background-color: #f5f5f5;
  border-top: 1px solid transparent;
  border-bottom-right-radius: 3px;
  border-bottom-left-radius: 3px;
}

.btn-info {
  color: #ffffff;
  background-color: #75caeb;
  border-color: #5fc1e8;
}
	")),
    shiny::tags$script(src = "shared/datatables/js/jquery.dataTables.min.js"), 
    shiny::tags$script(src = "shared/datatables/js/DT_bootstrap.js"),
    shiny::tags$script(src = "js/DTbinding.js"))),
    div(id = outputId, class = "shiny-datatable-output selectable"))
}
