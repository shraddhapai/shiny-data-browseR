suppressMessages(require(shiny))
suppressMessages(require(shinyBS))
suppressMessages(require(Gviz))
suppressMessages(require(Rsamtools))
suppressMessages(require(RColorBrewer))
suppressMessages(require(GenomicFeatures))

progressImage <- "images/spinner.gif"

#############################
# UI CONTROL OBJECTS . SEE BELOW for UI layout 
default.range <- c(10483022,10825218)

# #############################################################
# WELL PANELS (Plot settings and genome annotation)
# #############################################################
well_distStartPanel <- wellPanel(
style="border:1px inset;border-color:#458cc3;background-color:#ffffff",
  	HTML('<div class="cutefont">Genome Location</div><span style="font-style:italic" class="text-info">Manual plot update required</span><p>'),
    fluidRow(
        column(4,uiOutput("mychrom")),
        column(4,numericInput("xlim1", "x-range, from:", default.range[1])),
        column(4, numericInput("xlim2", "x-range, to:", default.range[2]))
    ),
    uiOutput("csize"),
    uiOutput("coordSize"),
    HTML('<p>&nbsp;</p><div class="cutefont"><b>Binning and Smoothing</b></div>'),
    fluidRow(
    	column(4,numericInput("nbin1", "Number of Bins:", 1000)),
    	column(4,uiOutput("bw2"))
    ),
    uiOutput("bwMess2")
)
 cat("\tWell: genomic location, bins, smoothing\n")
 
# well panel: choose data view type - individual, groups etc.,
well_chooseGroups <- wellPanel(
style="border:1px inset;border-color:#458cc3;background-color:#ffffff",
  HTML('<div class="cutefont">Data and Grouping</div><span class="text-info" style="font-style:italic">Manual plot update recommended</span><p>'),
  fluidRow(column(6, uiOutput("o_groupBy"))),
   fluidRow(column(6,selectInput("whichMetric", "Baseline trends:",
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

# well panel
well_yaxis <- wellPanel(   
style="border:1px inset;border-color:#458cc3;background-color:#ffffff",
  HTML('<div class="cutefont">Plot Options</div><span style="font-style:italic" class="text-info">Automatic refresh</span><p>'),
    fluidRow(
    	column(6,selectInput("plotType","Plot type", c("Mean+CI"="a_confint", "Points+Lines"="b", "Points"="p","Lines"="l","Smooth"="smooth"),"l"))	
    ), 
    fluidRow(
    	column(6,uiOutput("o_colorBy")),
    	column(6, selectInput("oCol", 
                              "Colour Scheme:", 
                              rownames(brewer.pal.info),"Spectral")
        )  
    ),
  	selectInput("whichYlim","Y-axis:", 
                  c("Default Bounds"="def","Custom Bounds"="cust"),"def"),
  conditionalPanel(
                   condition = "input.whichYlim == 'cust'",
                   div(class='row-fluid',uiOutput("customYlim"))
  ),
  checkboxInput("legd", "Show Legend", TRUE)
  )
cat("\tWell: Plot options\n")



# #############################################################
# shiny UI LAYOUT
# #############################################################
suppressWarnings(shinyUI(fluidPage(
  # fixed header bar
  HTML('<div class="header header-fixed">'),
  # main header showing title
  HTML('<div class="pagebar-col">'),
  fluidRow(
  	column(11, 
        HTML('<a href="https://github.com/shraddhapai/shiny-data-browseR"><img src="images/title.png", width=380></a>'),
        HTML('<span style="font-style:normal;font-size:14;vertical-align:-10px;align:left;font-weight:normal">The statistical browser for population ( and other ) genomics</span>')
    ),
    column(1,
  		conditionalPanel("input.getData>0",
  		shiny::tags$button(id="loadPlot2", type="button",
    	class="btn action-button btn-success", 
      HTML('<img src="images/refresh.png" style="width:20px;height:20px;align:right">'))
    ))  
  ),
  fluidRow(
   column(9,
          HTML('Version: <span style="color:#ffd357">EDB beta 1.1</span>')
   ),
   column(3,
          HTML('<span style="font-weight:800;font-size:16px;color:#ffd357">New!  <a href="https://www.youtube.com/watch?v=sv68ftm71R4" style="text-decoration:underline;color:#ffffff">"Quick Tour" screencast!</a>') 
   )),
  HTML('</div>'),
  HTML('</div>'),
  bsAlert('welcome_msg'),
    
    # #############################################################
    # COLLAPSIBLE PANELS
    # #############################################################
    bsCollapse(multiple=TRUE, open="col_activate",id="main_collapse",

       bsCollapsePanel(title="Choose dataset",
    	HTML('<div style="height:200px">'),
  		HTML('<h4 class="text-primary" style="margin-bottom:0px;">Select a dataset. Then click "Make active dataset".</h4>&nbsp;<br>'),
  		fluidRow(
    	column(7,uiOutput("pickData")),
    	column(2, shiny::tags$button(id="getData",
                                     type="button", class="btn action-button btn-lg btn-info", 
                                     HTML("Make active dataset")))
    	),
    	HTML("</div>"),
    	id="col_activate", value="activate"
    	),

    	bsCollapsePanel("Plot view",
    		conditionalPanel("input.getData>0",
               HTML('<div class="pagebar-col2">'),
               uiOutput("dataname"),
               fluidRow(
                  uiOutput("data_desc"),
                  column(1,
                            shiny::tags$button(id="loadPlot", type="button",
                            class="btn action-button btn-success", 
                            HTML('<img src="images/refresh.png" "style:width=10px;height:10px">')
                        )
               )),
               HTML("</div>"),
               bsAlert(inputId="plot_statusMsg"),
                 bsProgressBar("load_pBar", value=0,
                               visible=FALSE, color="standard",
                               striped=TRUE, animate=FALSE),
    			uiOutput("plot_welcome"),
    			 plotOutput("scatplot",height="400px")
    			),
    	id="col_plot",value="outplot"
    	),

    	bsCollapsePanel("Settings",
    	  	conditionalPanel("input.getData>0",
    	    HTML('<h4 class="text-primary" style="margin-bottom:0px;">Customize plot settings </h4>&nbsp;<br>'),
    		fluidRow(
    	  		column(4,well_distStartPanel),
    			column(4,well_chooseGroups),
    	 	   	column(4,well_yaxis)
    		)),
    	id="col_settings",value="settings"
    	), 

    	bsCollapsePanel("Sample selector",
  		 	conditionalPanel("input.getData>0",
    	HTML('<h4 class="text-primary" style="margin-bottom:0px;">Select samples for inclusion in analysis by multi-selecting rows below. </h4>&nbsp;<br>'),
  				HTML('<span class="text-info" style="font-style:italic">Manual plot update required</span><p>&nbsp;</p>'),
    			uiOutput("o_sampleCount"), 
    			selDataTableOutput("o_sampleTable")
    		),
    	id="col_sampleSel", value="sampleSel"
    	),

    	bsCollapsePanel("Genome Annotation",
    		conditionalPanel("input.getData>0",
    		HTML('<h4 class="text-primary" style="margin-bottom:0px;">Include genomic annotation as &quot;tracks&quot; below the main data.</h4><span style="color:#ff0000;font-style:italic">Note: Tracks showing gene models can add >1min to refresh time, depending on view range.<br>It is recommended that these be turned on at the end after plot settings have been finalized.</span>&nbsp;<br>&nbsp;<br>'),
  				HTML('<span class="text-info" style="font-style:italic">Manual plot update required</span>'),
   				wellPanel(uiOutput("o_getAnnot"))
    	),
    	id="col_genomeAnnot", value="genomeAnnot"
    	)
    ), # end bsCollapse

    # #############################################################
    # UI TOOLTIPS
    # #############################################################
    bsTooltip("dataset", 
              'Select a dataset to analyze. Then click "Make active dataset".',
              "top"),
    # Settings: Genomic Location
    bsTooltip("mychrom",
              "Select genomic coordinates. Only available sequences are shown.", 
              "top"),
    bsTooltip("nbin1",
              "Binning the data - how many bins for viewable range?", 
              "top"),
    bsTooltip("bw2", 
              "Bandwidth for Gaussian smooth. 50% of the data lies within a quarter of this value",
              "right"),
    # Settings: Data and grouping
    bsTooltip("whichMetric",
              "Divide sample/group trends by a baseline.", 
              "right"),
    # Settings: Plot Options
    bsTooltip("plotType", 
              "Plot trendlines using points, lines or with errorbars",
              "right"),
    bsTooltip("colorBy", 
              "Color trendlines by categorical variable. Only works when 'Group samples by' is set to '(none)'",
              "top"),
    bsTooltip("oCol",
              "Pick color palette for colouring trendlines. See colorbrewer2.org for descriptions",
              "top"),
    bsTooltip("whichYlim", 
              "Select if y-axis should have default or custom limits",
              "right"),
    fluidRow(
     column(7,
        HTML('<h4>For code, technical documentation, and user manual, visit <a href="https://github.com/shraddhapai/shiny-data-browseR">EDB @ github</a></h4>')
    ), column(1,HTML("")
    ), column(4,
    HTML('<h5 style="text-align:right">Report bugs to: Shraddha | dot | Pai | at | camh | dot | ca</h5>')
    )),
    HTML('<div class="panel-footer">Copyright &copy; 2014. Shraddha Pai, Centre for Addiction and Mental Health. This software is distributed with the GPLv3 license.<br>This page is best viewed at 1280 x 1024 or better, and has been tested in Firefox 31.0 and Chrome 37.0 </div>')
    )))
