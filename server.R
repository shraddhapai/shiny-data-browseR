suppressMessages(require(shiny))
suppressMessages(require(GenomicRanges))
suppressMessages(require(RColorBrewer))
suppressMessages(require(rtracklayer))
suppressMessages(require(doMC)); 

registerDoMC(5) #TODO make this a param in config file

source("process_samps.R")
source("plotters.R")

# defaults
verbose <- T #### set to T for debugging
if (verbose) cat("************\nDebug mode\n*****************\n")
# increase max file size for reference chromosomes to 10MB
options(shiny.maxRequestSize=10*1024^2)
options(scipen=10)

configSet <- listConfig() # load config for all datasets 

# Define server logic required to plot variables
shinyServer(
function(input, output, session){
		cat("In shinyServer\n")
		updateCollapse(session, id = "main_collapse",  
                       open = "col_plot", close = NULL)

	# which button was last pressed?
	# this is needed to separate data load from plot plot.
	values <- reactiveValues()
	values$lastAction <- NULL
	observe({if ( input$getData!=0 ) { values$lastAction <- "data"}})
	observe({if ( input$loadPlot!=0 || input$loadPlot2!=0 ) {values$lastAction <- "plot"}})

	# tier 0 : shows up when page is loaded.
	output$pickData <- renderUI({selectInput("dataset", "", names(configSet), width="750px")	})
	# tier 1 : happens when 'make active dataset' button is clicked.
    refreshConfig <- reactive({
		plot_alertTxt <- paste(
			"<div style='font-size:20px;font-weight:800;margin-bottom:10px'>",
            'Click <img src="images/refresh.png" style="width:20px;background-color:#00cc00;margin:5px; padding:3px"> to see a default plot.<br>',
            "Or scroll down to customize first.<br>",
            "</div>",
			sep="")

   		if (input$getData == 0) return(NULL) # only depends on first button
		createAlert(session, inputId="plot_statusMsg", 
                    alertId="alert_statusMsg", message=plot_alertTxt, 
                    type="info",dismiss=FALSE,append=FALSE)
		updateCollapse(session, id="main_collapse", open="col_settings", 
                       close="col_activate")
		
        if (verbose) cat("* Refreshing config")
		settings <- configSet[[input$dataset]]
		allDat <- settings$allDat
		chromsize <- settings$chromsize
		groupKey <- settings$groupKey
		configParams <- settings$configParams

		isolate({return(settings)})
   })

  output$o_groupBy <- renderUI({ 
  if (input$getData == 0) return(NULL)
  	settings <- isolate({refreshConfig()}); if (is.null(settings)) return(NULL) 
	groupNames <- settings$configParams$groupCols
	selectInput("groupBy", "Group samples by:", c(groupNames,"(none)"),
		selected=settings$configParams$defaultGroup)
})
  #if (input$getData > 0) {
#	addTooltip(session, id="groupBy","Generate mean trendlines by categorical variable. (none) for see individual samples as-is","right")
 # }
  if (verbose) cat("\tGot by groupBy\n")

  output$dataname <- renderUI({
	blank <- HTML(paste('<div style="height:50px;font-style:italic;color:#c0e4ff;margin-left:10px;margin-top:10px">',
                        '<span style="font-size:20px;font-style:normal;color:#ffd357;font-weight:800">Welcome!</span>',
                        '<br>The view below is composed of 5 panels. Each panel can be collapsed or expanded by clicking ',
                        'on the respective titles.<br>Begin data exploration by selecting a dataset in the first panel, ',
                        '"Select dataset".</div>',sep=""));
   if (input$getData == 0) return(NULL)
  	settings <- isolate({refreshConfig()}); if (is.null(settings)) return(blank) 
	fluidRow(
		column(9,
			HTML(paste('<div style="color:#ffffff;font-size:18px;margin-top:10px;margin-left:10px">Active dataset:',
			sprintf('<span style="color:#ffd357; font-weight:600">%s</span>', settings$configParams[["name"]])),
			sprintf(': build <span style="color:#ffd357;font-weight:600">%s</span></div>', 
			settings$configParams[["genomeName"]]),
			sep="")
	), 	column(3,HTML(sprintf('<em style="margin-top:10px;color:#ffd357">%i samples available</em>', 
		nrow(settings$allDat))))
	)})


  output$data_desc <- renderUI({
	blank <- HTML('<div height:50px">&nbsp;</div>');
   if (input$getData == 0) return(blank)
  	settings <- isolate({refreshConfig()}); if (is.null(settings)) return(blank) 
	fluidRow(
		column(8,HTML(sprintf('<div style="margin-left:10px;margin-top:10px"><span style="#ffffff;font-weight:800">Description:  </span><span style="color:#c0e4ff;font-weight:400">%s</span></div>',settings$configParams[["description"]]))),
		column(1,HTML("")),
		column(2,HTML(sprintf('<span>Platform/assay:<br><span style="color:#c0e4ff">%s</span></span>', 
		settings$configParams[['platformName']])))
	)
	})

	output$o_sampleCount <- renderUI({
  	settings <- isolate({refreshConfig()}); if (is.null(settings)) return(NULL) 
	if (input$getData==0) return(NULL)
	if (is.null(input$o_sampleTable)) {
		x <- nrow(settings$allDat) 
	} else if (length(input$o_sampleTable)==1) {
		if (input$o_sampleTable==-1) {
		x <- 0 #nrow(settings$allDat) 
		}
	}  else {
		tmp <- matrix(input$o_sampleTable,byrow=T,ncol=ncol(settings$allDat)-1); 
		x <- nrow(tmp)
	}
	HTML(sprintf("<i>%i samples selected</i>\n", x))
	})
  if (verbose) cat("\tGot by sampleCount\n")
  
  ### region help text
  output$mychrom <- renderUI({
   if (input$getData == 0) return(NULL)
  	settings <- isolate({refreshConfig()}); if (is.null(settings)) return(NULL) 
    
   source(sprintf("data_types/%s.R", settings$configParams$datatype))
    myfiles <- settings$allDat; sq <- getSeqinfo(myfiles$bigDataURL[1])
    selectInput("chrom", "Sequence:", sq)
  })
  
  output$csize <- renderUI({ #chrom size message
   if (input$getData == 0) return(NULL)
  	settings <- isolate({refreshConfig()}); if (is.null(settings)) return(NULL) 
    if(!is.null(input$chrom)){
    thischromsize <- settings$chromsize[which(settings$chromsize$chrom == input$chrom),2]
    HTML(sprintf("<b>Chrom max: %s bp</b>", format(thischromsize,big.mark=",")))
    }
  })
  
  output$coordSize <- renderUI({ #region width message
    r1 <- input$xlim1; r2 <- input$xlim2
    HTML(sprintf("<b>Region Width: %6.1f kb, %3.1f Mb</b>", (r2-r1)/1e3,(r2-r1)/1e6))
  })
  
  #### choosing data:
  output$customYlim <- renderUI({
	myout <- refreshData(); 
	if (is.null(myout))  return(NULL) 
	else {
		outdat <- myout[["outdat"]]
		tmp <- na.omit(as.numeric(as.matrix(outdat[,-(1:3)]))); 
		values <- c(quantile(tmp,c(0.005,0.995)),min(tmp),max(tmp))
	}
	print(values)
	sliderInput("customYlim", "Custom y-range", value=values[1:2],min=values[3],max=values[4])
  })

  output$o_colorBy <- renderUI({
  if (input$getData == 0) return(NULL)
  settings <- isolate({refreshConfig()}); if (is.null(settings)) return(NULL) 
  if (is.null(input$groupBy)) g <- settings$configParams$defaultGroup else g <- input$groupBy
	groupNames <- settings$configParams$groupCols
	selectInput("colorBy", "Color by:", c(groupNames,"(none)"),selected=g)
  })
  
  output$out_baseline <- renderUI({
   if (input$getData == 0) return(NULL)
  	settings <- isolate({refreshConfig()}); if (is.null(settings)) return(NULL) 
	if (is.null(input$groupBy)) g <- settings$configParams$defaultGroup else g <- input$groupBy
	if (g=="(none)") selVal <- NA else selVal <- settings$groupKey[[g]]
  	selectInput("whichBaseline","Baseline by", selVal)
  })
  
  
  #### for smoothing:
  output$bw2 <- renderUI({
      if(!is.na(input$xlim1) & !is.na(input$xlim2)){
        r1 <- input$xlim1
        r2 <- input$xlim2
        default.param <- (diff(c(r1,r2)))*0.02 # for default 5%
        numericInput("param2", "Smooth bw (bp):",
                     min=0, value=as.integer(default.param))
      }
  })
  
  output$bwMess2 <- renderUI({
    if(!is.na(input$xlim1) & !is.na(input$xlim2) & !is.null(input$param2)){
    rangeSize <- (input$xlim2-input$xlim1)
    bwSize <- input$param2
    bwProp <- (bwSize/rangeSize)*100
	binSize <- rangeSize/input$nbin1
    HTML(sprintf("<b>Bin size = %6.1f kb or %3.1f Mb<br>Bandwidth = %2.2f%% of win</b>", 
		binSize/1e3, binSize/1e6,bwProp))
    }
  })


  # sample table
  output$o_sampleTable <- renderDataTable({
  	if (input$getData==0) return(NULL)
  	settings <- isolate({refreshConfig()}); if (is.null(settings)) return(NULL)
	df <- settings$allDat; df <- df[,-which(colnames(df)=="bigDataURL")]
	return(df)
  }) #,options=list(bSortClasses=TRUE))

  # annotation view
  output$o_getAnnot <- renderUI({
  if (input$getData==0) return(NULL)
  	settings <- isolate({refreshConfig()}); if (is.null(settings)) return(NULL)
	anno <- read.delim(settings$configParams$annoConfig,sep="\t",header=T,as.is=T)
	checkboxGroupInput("anno", "The following annotation tracks are available for the current genome build:", choices=anno$name, selected=NULL)
  })

    # ######################################################################
    # TIER 2 : happens when 'Compute Plots' button is clicked.
    # ######################################################################
	refreshData  <- reactive({
    
        if(input$loadPlot == 0 && input$loadPlot2 == 0) return(NULL)
        return(isolate({ # everything in this function waits for actionButton() 
                         # to be selected before executing
    		if (verbose) cat("* In refreshData()\n")
    
    		settings <- refreshConfig()
            cat("\t got config\n")
    
    		# re-create data matrix from selectableDataTable input
    		myfiles <- settings$allDat
    	
    		if (is.null(input$o_sampleTable)) {
    			idx <- 1:nrow(myfiles)
    		} else if (length(input$o_sampleTable)==1 && input$o_sampleTable==-1) {
    			idx <- 1:nrow(myfiles)
    		} else {
    			tmp <- matrix(input$o_sampleTable,byrow=T,ncol=ncol(myfiles)-1); 
    			colnames(tmp) <- colnames(myfiles)[-which(colnames(myfiles)=="bigDataURL")]
    			prefilter_samps <- tmp[,"sampleName"]
    			idx <- which(myfiles$sampleName %in% prefilter_samps)
    		}
    		if (!any(idx)) return("Please select at least 1+ sample using the Sample Selector");
    
    		if(length(idx)==1) {tmp <- colnames(myfiles);  myfiles <- as.data.frame(myfiles[idx,]); colnames(myfiles) <- tmp;}
    		else { myfiles <- myfiles[idx,]}
    		cat(sprintf("After sample filtering, have %i samples\n", length(idx)))
    
            if(length(myfiles$sampleName) != unique(length(myfiles$sampleName))){ cat("There are repeat samples!\n"); browser()}
    
           cat("* Read files, bin data, create matrix\n") 
    	   selRange <- GRanges(input$chrom, IRanges(input$xlim1,input$xlim2))
    	    
    	   createAlert(session, inputId="plot_statusMsg", alertId="alert_statusMsg", 
                       message="Fetching data", type="warning",dismiss=FALSE,append=FALSE)
    	   print(system.time(dat <- fetchData(session, myfiles=myfiles, selRange=selRange,
                                            numBins=input$nbin1,
                                            datatype=settings$configParams$datatype,
                                            datatypeParams=settings$datatypeParams,
                                            verbose=verbose)))
    	   updateProgressBar(session, "load_pBar", visible=FALSE)
    	   bed <- dat$coords; alldat <- dat$values; rm(dat)
    
    		baselineTxt <- ""
    		if (input$groupBy != "(none)") {
    			if (input$whichMetric != "normal") {
    	   			createAlert(session, inputId="plot_statusMsg", alertId="alert_statusMsg",
                                   message="Baselining samples", type="warning",dismiss=FALSE,
                                   append=FALSE)
            		mygroups <- settings$groupKey[[input$groupBy]]
    				mygroups <- intersect(mygroups, myfiles[,input$groupBy])
    				
    				grp.summ <- computeAverages(myfiles,mygroups,bed,alldat,F,input$groupBy)
    				cat("\taverage computed\n")
    				tmp <- baselineSamps(alldat, grp.summ, mygroups, 
    						input$whichMetric,input$whichBaseline, input$logMe)
    				alldat <- tmp$alldat; 
    				logTxt <- "unlogged"; if (input$logMe) logTxt <- "log2"
    				baselineTxt <- sprintf(": %s:%s (%s)", input$whichMetric, input$whichBaseline, 
                                           logTxt)
    			}
    		}
    		outdat <- cbind(bed,alldat)
    
            myOut <- list(myfiles=myfiles, outdat=outdat,plotTxt=sprintf("%s %s", settings$configParams$name, baselineTxt))

    		cat("\tReturning output\n")
            
            return(myOut)
    		### 1) myfiles: phenotype table
    		### 2) outdat: data.frame with 3+N columns, where N is number of samples. First three columns are named chrom,start,end.
    		### 3) plotTxt: character for plot title
    		### 4) baselineTxt: title suffix indicating if samples have been baselined. -- OBSOLETE?
          })) # end of return(isolate()) 
  })

  output$scatplot <- renderPlot({ 
    while (sink.number() >=1) { sink(NULL)} # shinyBS seems to open log files. Close these.
    
    cat("In renderPlot\n")

    # cases where plot doesn't have all information to load the plot
    if (is.null(values$lastAction)) return(NULL)
    if (values$lastAction=="data") return(NULL)
    if (input$loadPlot==0 & input$loadPlot2 == 0 ) return(NULL)
    
    updateCollapse(session, id = "main_collapse", multiple=TRUE,
                   close="col_activate", open="col_plot")

    myOut <- isolate({refreshData()});if(is.null(myOut)) return(NULL)

    cat("past refreshData\n")
    settings <- isolate({refreshConfig()})  
    cat("past refreshConfig")
    
    if(class(myOut) == "character") return(NULL) # print no graph
    else{
        cat("* Render plot\n")
        myfiles <- myOut[["myfiles"]]
        outdat <- myOut[["outdat"]]
        settings$groupKey[["sampleName"]] <- myfiles$sampleName
    
        # update available groups based on selected samples.
        groupKey <- settings$groupKey
        for (k in names(groupKey)) {
            x <- groupKey[[k]]
            #TODO: define behaviour if no items match here
            groupKey[[k]] <- x[which(x %in% unique(myfiles[,k]))]
        }

        #TODO what if groupBy has no entries left?
        if (input$groupBy !="(none)" && is.null(groupKey[[input$groupBy]])){
            cat("you haven't decided what to do here, SP")
        }
    
        # defer plot refresh if these params change
        isolate({Smoother2 <- input$param2})
        isolate({selAnno <- input$anno})
        isolate({xlim <- c(input$xlim1, input$xlim2)})
        closeAlert(session,alertId="alert_statusMsg")
        
        mkScat(session, statusId="plot_statusMsg", myfiles=myfiles, 
            outdat = outdat,configParams=settings$configParams,
            xlim = xlim,
            groupKey=groupKey, groupBy=input$groupBy, 
            colorBy=input$colorBy, oCol=input$oCol, 				# color
            plotViewType=input$plotType, plotType="smoo2", 		# plot type
            param2=Smoother2, 			# errorbar related
            whichYlim=input$whichYlim, customYlim=input$customYlim,	# ylim
            legd=input$legd,											# legend 
            plotTxt=myOut$plotTxt,selAnno=selAnno,
            verbose=TRUE
        )
    
        updateCollapse(session, id="main_collapse", 
        open=c("col_plot", "col_settings"))
    }
  })

}) # end of shinyServer function


