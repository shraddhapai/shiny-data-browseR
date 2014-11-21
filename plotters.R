## Functions which produce Plots
mkScat <- function
(
session, 
statusId, 
myfiles, # the config file
outdat, # the main output from reactive function
configParams, ##<<(list) key-value pairs from config.txt
xlim,
groupKey,	##<<(list) key: group category, value: order of group names
groupBy,	##<<(character) column to group by. Must be a key in "groupKey" object
colorBy,	##<<(character) variable to colour by. Must be a key in "groupKey" object
oCol, 		##<<(character) group colour scheme
plotViewType, 	##<<(character) plot type for Data Track. See Gviz for options
plotType, 	##<<(character) currently hard-coded to "smoo2"
param2, # smoothing bandwidth
whichYlim, # whether default or custom Y axis was chosen
customYlim,
legd, # T/F whether legend should be shown
plotTxt,
selAnno,
verbose=T # for debugging
){

if (verbose) "In plotters::mkScat\n"

tryCatch({

if("chrom" %in% colnames(outdat)){
	setchrom <- unique(as.character(outdat$chrom))
	outdat <- outdat[,which(!(colnames(outdat) == "chrom"))]
}
outdat <- outdat[order(outdat$start),]


if (verbose) cat("\toutdat successfully reordered\n")

# defining variables
if (groupBy=="(none)") {
	mygroups <- groupKey[["sampleName"]]
	if (colorBy=="(none)") { 
            colgroups <- groupKey[["sampleName"]]; 
            colorBy <- "sampleName" 
    }
	else colgroups <- groupKey[[colorBy]]
} else {
	mygroups <- groupKey[[groupBy]]
	colgroups <- mygroups
}

# x/y limits, smoothing param
if(whichYlim == "def"){
	tmp_mat <- as.matrix(outdat[,-c(1,2,grep("min",colnames(outdat)), grep("max",colnames(outdat)))])
	tmp <- na.omit(as.numeric(tmp_mat)); myYlim <- quantile(tmp,c(0.005,0.99)); 
} else myYlim <- customYlim

bw <- param2
xvals <- outdat$start

if (verbose) cat("\tset up xylim successfully\n")

# colours
colBase <- suppressWarnings(brewer.pal(length(colgroups)+1,oCol)[-1])
if (length(colBase) > length(colgroups)) {
        colBase <- colBase[1:length(colgroups)]
}
colist <- colBase
coldiff <- length(colgroups)-length(colist)
while (coldiff > 0) { 
	colist <- c(colist, colist[1:min(coldiff, length(colBase))]) # wrap around when n > num colours
	coldiff <- length(colgroups)-length(colist)
	}
rm(colBase)
names(colist) <- colgroups # want colours to stay the same regardless of which groups being plotted
#pdf("~/test.pdf"); barplot(1:length(colist), col=colist, border=NULL); dev.off()
#print(colist)

cat("\t* Setting up annotation tracks\n")
gtrack <- GenomeAxisTrack(cex=1.5)
itrack <- getIdeo(configParams$ideoFile, configParams$genomeName, setchrom)
plot_GR <- GRanges(seqnames=setchrom, IRanges(outdat$start, outdat$start))

# ####################################
# View mode 1: Group-wise plot
if(!groupBy %in% "(none)"){ 
	if (is.null(mygroups)) return(NULL)
	if (plotViewType=="l") plotViewType <- "a" # when grouping, we want to show the mean
	
    if (ncol(outdat) == 3) { # single sample selected 
        tmp_name    <- colnames(outdat)[3]
        outdat      <- outdat[,-(1:2)]
        outdat      <- as.matrix(outdat)
        colnames(outdat) <- tmp_name
    } else {
	    outdat <- outdat[,-(1:2)]
    }
	
	# set colour and fill
	filldat <- colist
	if ("a_confint" %in% plotViewType){
		source("makeColorLighter.R")
		plotViewType	<- c("a","confint")
		filldat <- sapply(colist, makeColorLighter)
	} 
	
	# smoothing if necessary
	vals <- outdat
	if (plotType=="smoo2") {
		vals <- smooMat(outdat,bw,xvals,1:ncol(outdat))
	}
	mcols(plot_GR) <- vals;
    #head(plot_GR)
    #print(summary(vals))
	
	# create data track
	cat("\t* Created data track\n")
	g <- groupKey[[groupBy]]; g <- g[which(g %in% mygroups)]
	cat(sprintf("Plot view type=%s\n", plotViewType))

    g_samples   <- factor(myfiles[,groupBy],levels=g);
    g_count     <- as.integer(table(g_samples)) 
    
    # convert groupname from "myGroup" to "myGroup (5)" where 5 is the number of samples in that group.
    g_newGroupNames  <- paste(g, paste(paste(rep("(",length(g_count)), g_count,sep=""), rep(")",length(g_count)),sep=""),sep=" ")
    g_newSampGroups <- myfiles[,groupBy]
    for (i in 1:length(g)) { 
            g_newSampGroups[which(g_newSampGroups == g[i])] <- g_newGroupNames[i]
    }

	dTrack <- DataTrack(plot_GR,groups=factor(g_newSampGroups,levels=g_newGroupNames),levels=g_newGroupNames,
		col=colist,fill=filldat, type=plotViewType,name=configParams[["ylabel"]])
	
	# ######################################################
	# View mode 2: Individual samples
} else{ 
	# reorder based on group priority
  
	smpsGrps <- getSampleOrder(myfiles,groupKey)
	midx <- match(myfiles$sampleName, colnames(outdat))
	if (all.equal(colnames(outdat)[midx],myfiles$sampleName)!=T) { 
            cat("sample names don't match output matrix column names. Check?")
    }
	outdat <- outdat[,midx]
	if (nrow(myfiles)==1) { outdat <- as.matrix(outdat); colnames(outdat) <- myfiles$sampleName }

	# assign sample colours
	colist <- colist[match(smpsGrps[,colorBy], names(colist))]
	
	if (plotType == "smoo2") {
            outdat <- smooMat(outdat,param2, xvals,1:ncol(outdat))
    }
	mcols(plot_GR) <- outdat
	dTrack <- DataTrack(plot_GR, groups=factor(smpsGrps[,"sampleName"],levels=smpsGrps[,"sampleName"]),
		aggregateGroups=F,type=plotViewType,
		col=colist,name=configParams[["ylabel"]])
} 


# ######################################################
# Track setup and plot

cat("setting display track\n")
displayPars(dTrack)  <- list(
	cex.title=1.5,background.title='cyan4', 
	#baseline=0, col.baseline='gray20',lty.baseline=1,
    legend=legd,cex.axis=1.5,lwd=2,
	col.grid='gray50',lty.grid=3,cex=1,pch=20, cex.legend=1.3,
	grid=TRUE,ylim=myYlim
	)

# Getting annotation
createAlert(session, inputId=statusId, alertId="alert_statusMsg", 
            message="Fetching annotation", type="warning",dismiss=FALSE,append=FALSE)
cat("* About to get anno\n")
sizes <- c(0.05,0.18,1.0)
finalTracks <- list(itrack, gtrack, dTrack)
if (!is.null(selAnno)) {
	cat("\t* Collecting annotation tracks\n")
	annoTracks <- getAnnoTracks(configParams$annoConfig, selAnno,
                                GRanges(setchrom,IRanges(xlim[1],xlim[2])),configParams$genomeName)
    cat("GOt out of getAnnoTracks\n")
	finalTracks <- c( finalTracks, annoTracks[["tracks"]])
    cat("Made finalTracks\n")
	sizes <- c(sizes, annoTracks$sizes)
}

# Create plot title string
groupStr <- ""
nSamp <- ncol(mcols(plot_GR))
if (groupBy == "(none)") {
    groupStr <- sprintf("%i Samples, UNGROUPED", nSamp)
} else {
    groupStr <- sprintf("%i Samples, AVERAGED BY %s", nSamp, groupBy)
}

if (diff(xlim) > 2e6) {
    xlim_str <- sprintf("%1.1f - %1.1f Mb", xlim[1]/1e6, xlim[2]/1e6)
} else {
    xlim_str <- sprintf("%s - %s bp", 
                        prettyNum(round(xlim[1]/1e3),big.mark=","), 
                        prettyNum(round(xlim[2]/1e3),big.mark=",")
                        )
}

ttl <- sprintf("%s: %s: %s\n%s", plotTxt,setchrom, xlim_str, groupStr)
createAlert(session, inputId=statusId, alertId="alert_statusMsg", 
            message="Rendering tracks", type="warning",dismiss=FALSE,
            append=FALSE)

eleft <- max(10,0.1 * diff(xlim))

cat("about to call plotTracks\n")
t0 <- system.time(
	plotTracks(finalTracks, sizes=sizes,
		min.width=1, extend.left=eleft, lwd=2,
		fontface.main=2,cex.axis=1.5,
		chromosome=setchrom,from=xlim[1], to=xlim[2],
		main=ttl,cex.main=1.2
	)	
)
if (verbose) {cat("Plotting time:\n"); print(t0) }
}, error=function(ex) {
	cat("Error in plotters::mkScat:\n")
	print(ex)
}, finally={
})
cat("--------\n")
closeAlert(session,alertId="alert_statusMsg")
}


#########################################################
# Functions to create Gviz annotation tracks
#########################################################

getAnnoTracks <- function#
(anno, selAnno,selGR,genomeName) {
dat <- read.delim(anno,sep="\t",header=T,as.is=T)
dat <- dat[which(dat$name %in% selAnno),] 
cat("In getAnnoTracks\n")

registerDoMC(6)
print(system.time(trackList <- foreach (k=1:nrow(dat)) %dopar% {
	curr <- dat[k,]
	gr <- switch(tolower(curr$format),
		tabix= {
			tbx <- scanTabix(TabixFile(curr$bigDataURL), param=selGR)	
			recs <- strsplit(tbx[[1]],"\t")
			mat <- matrix(unlist(recs), ncol=length(recs[[1]]), byrow=T)
			rm(recs)
			GRanges(mat[,1],IRanges(as.integer(mat[,2]), as.integer(mat[,3])))
		},
		bigwig= {
			import.bw(curr$bigDataURL,asRangedData=F,which=selGR)
		},
		txdb={
			src <- curr$bigDataURL
			if (any(grep("^TxDb", src))) { require(src); txdb <- eval(parse(text=src)) 
			} else { txdb <- loadDb(src)
			}
			seqlevels(txdb,force=TRUE) <- seqlevels(selGR)
		},
		stop(sprintf("Invalid 'format' for track: %s (format given = %s)", curr$name))
	)	
	tr <- switch(curr$trackType, 
		AnnotationTrack={
			x <- AnnotationTrack(gr, name=curr$name, stacking=curr$defaultView,
				genome=genomeName, chrom=seqlevels(selGR))
			x
		}, 
		GeneRegionTrack={
			if ((end(selGR)-start(selGR)) < 4e5) { 
				viewOpt <- curr$defaultView; showDetail <- TRUE
			} else {
				viewOpt <- "dense"; showDetail <- FALSE
			}
			x <- GeneRegionTrack(txdb, chromosome=seqlevels(selGR), genome=genomeName,
				stacking=viewOpt, name=curr$name, showId=showDetail, geneSymbol=showDetail,
				collapseTranscript=!showDetail,
				shape="smallArrow")
			x	
		},
		stop(sprintf("Unsupported trackType: %s", curr$trackType))
		)
		displayPars(tr) <- list(cex.title=1.2, col.title="white",background.title="darkred", rot.title=90,
			fill=curr$color, col=curr$color)
		return(tr)
}))

return(list(tracks=trackList,sizes=dat$sizes))
}



getIdeo <- function# return ideogram track
( 
inFile, ##<<(character) ideogram file from config.txt
genomeName, ##<<(character) genome name
chr
) {
    dat <- read.delim(inFile,sep="\t",header=T,as.is=T)
    
    return(IdeogramTrack(chr, genomeName, bands=dat))
}


