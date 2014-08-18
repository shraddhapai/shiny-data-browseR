suppressMessages(require(shiny))
suppressMessages(require(shinyBS))
suppressMessages(require(Gviz))
suppressMessages(require(Rsamtools))
suppressMessages(require(RColorBrewer))
suppressMessages(require(GenomicFeatures))

progressImage <- "images/spinner.gif"
tooltipStr <- '<script type="text/javascript" src="js/opentip-jquery.min.js"></script><link href="css/opentip.css" rel="stylesheet" type="text/css"><script type="text/javascript" src="js/library.js"></script>';

#############################
# UI CONTROL OBJECTS . SEE BELOW for UI layout 
# Defaults
#default.range <- c(93305559,93359905)
default.range <- c(0,5e5) + 20e6

# well - dist start panel
well_distStartPanel <- wellPanel(
style="border:1px inset;border-color:#458cc3;background-color:#ffffff",
  	HTML('<div class="cutefont">Genome Location</div><i>Manual button refresh required</i><p>'),
    fluidRow(
        column(4,uiOutput("mychrom")),
        column(4,numericInput("myrange1", "x-range, from:", default.range[1])),
        column(4, numericInput("myrange2", "x-range, to:", default.range[2]))
	),
	bsTooltip("myrange1", "Select range of x-axis view"),
    uiOutput("csize"),
    uiOutput("coordSize"),
	HTML('<p>&nbsp;</p><div class="cutefont"><b>Binning and Smoothing</b></div>'),
	fluidRow(
		column(4,numericInput("nbin1", "Number of Bins:", 1000)),
		column(4,uiOutput("bw2"))
	),
	bsTooltip("nbin1", "Select number of bins to group x-axis data in. Width of a bin shown below","right"),
	bsTooltip("param2","Select bandwidth for smoothing data. Fraction of window covered by bandwidth shown below.","right"),
	uiOutput("bwMess2")
)
 cat("\tWell: genomic location, bins, smoothing\n")
 
# well panel: choose data view type - individual, groups etc.,
well_chooseGroups <- wellPanel(
style="border:1px inset;border-color:#458cc3;background-color:#ffffff",
  HTML('<div class="cutefont">Data and Grouping</div><div style="color:red"><i>Manual button refresh recommended</i></div><p>'),
  fluidRow(column(6, uiOutput("o_groupBy"))),
  bsTooltip("groupBy", "Select sample property to average by","top"),
   fluidRow(column(6,selectInput("whichMetric", "Baselining option:", 
			c("Sample-wise data (no baselining)"="normal", 
			"Sample/Mean-of-all"="ratio_all", 
			"Sample/Mean-of-baseline"="ratio_baseline"),"normal")),
	conditionalPanel(
			condition = "input.whichMetric != 'normal'",
			column(6,checkboxInput("logMe","log-transform: log(1+x,base=2)", TRUE))
	)),
	conditionalPanel(
    		condition = "input.whichMetric == 'ratio_baseline'",
			fluidRow(column(8,uiOutput("out_baseline")))
		)
)
 cat("\tWell: Choose groups\n")


well_yaxis <- wellPanel(   
style="border:1px inset;border-color:#458cc3;background-color:#ffffff",
  HTML('<div class="cutefont">Plot Options</div><div style="color:red"><i>Automatic refresh</i></div><p>'),
    fluidRow(
		column(6,selectInput("plotType","Plot type", c("Mean+CI"="a_confint", "Points+Lines"="b", "Points"="p","Lines"="l","Smooth"="smooth"),"l"))	
	), 
	fluidRow(
		column(6,uiOutput("o_colorBy")),
		column(6, selectInput("oCol", "Colour Scheme:", rownames(brewer.pal.info),"Dark2"))  
	),
  	selectInput("whichYlim","Y-axis:", c("Default Bounds"="def","Custom Bounds"="cust"),"def"),
  conditionalPanel(
    condition = "input.whichYlim == 'cust'",
    div(class='row-fluid',uiOutput("customYlim"))
    ),
  	checkboxInput("legd", "Show Legend", TRUE)
  )
 cat("\tWell: Plot options\n")
 
# --------------------------------------------------------------------
# Main UI layout
# --------------------------------------------------------------------
   HTML_ButtonClickText <-   HTML(paste('<script type="text/javascript">',
       '$(document).ready(function() {',
          '$("#loadPlot").click(function() {',
            sprintf('$("#scatplot").html(\"<h4>Computing plot.<br>Do not change settings till plot has refreshed.</h4><img src=\\"%s\\">\");', progressImage),
          '});',
       '});',
      '</script>',
	  sep=""))

            # :'$("#errorMessage1").text("");',
suppressWarnings(shinyUI(fluidPage(#theme="bootstrap.css",
  HTML('<div class="header header-fixed pagebar-col">'),
  fluidRow(
  	column(9, HTML('<div>&nbsp;</div><span style="font-family:\'Homemade Apple\',cursive; font-size:24px">&nbsp;Epigenome Data Browse-R</span><div>&nbsp;</div>')),
	column(3,
  		conditionalPanel("input.getData>0",
  		shiny::tags$button(id="loadPlot", type="button",
		class="btn action-button btn-success btn-xlg", HTML("   Update Plot   "))
	))),
  HTML('<div class="cutefont,subtitle" style="font-size=200%">'),uiOutput("dataname"),HTML('</div>'),
  HTML('</div>'),
  HTML_ButtonClickText,
	bsCollapse(multiple=TRUE, open="col_activate",id="main_collapse",
	   bsCollapsePanel(title="Activate dataset",
		HTML('<div style="height:300px">'),
  		fluidRow(
  		column(3,HTML('<h4 class="text-primary" style="margin-bottom:0px;">Initial: Activate a dataset</h4>')),
		column(7,uiOutput("pickData")),
		column(2,
			shiny::tags$button(id="getData",
							   type="button", class="btn action-button btn-lg btn-info", HTML("Make active dataset")))
		),
		HTML("</div>"),
		id="col_activate", value="activate"
		),
			bsCollapsePanel("Plot",
				conditionalPanel("input.getData>0",
								 plotOutput("scatplot",height="300px")
				),
			id="col_plot",value="outplot"
			),
			bsCollapsePanel("Settings",
			  	conditionalPanel("input.getData>0",
				fluidRow(
			  		column(4,well_distStartPanel),
					column(4,well_chooseGroups),
			 	   	column(4,well_yaxis)
				)),
			id="col_settings",value="settings"
			), 
			bsCollapsePanel("Sample selector",
  		 		conditionalPanel("input.getData>0",
				HTML('<div style="background-color:#cccccc"><b>Sample Selector</b>'),
				uiOutput("o_sampleCount"), selDataTableOutput("o_sampleTable"),
				HTML("</div>")),
			id="col_sampleSel", value="sampleSel"
			),
			bsCollapsePanel("Genome Annotation",
				conditionalPanel("input.getData>0",
  				HTML('<span style="color:red"><i>Manual button refresh required</i></span>'),
   				wellPanel(fluidRow(uiOutput("o_getAnnot"))),
				HTML("</div>")),
			id="col_genomeAnnot", value="genomeAnnot"
			)
	), # end bsCollapse
			bsTooltip("plotType", "Select a plotType","right"),
			HTML('<div class="panel-footer">Copyright &copy; 2014. Shraddha Pai. This software is distributed with the GPLv3 license. </div>')
	)))
