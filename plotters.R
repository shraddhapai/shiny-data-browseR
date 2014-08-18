## Functions which produce Plots
mkScat <- function
(
myfiles, # the config file
outdat, # the main output from reactive function
configParams, ##<<(list) key-value pairs from config.txt
groupKey,	##<<(list) key: group category, value: order of group names
groupBy,	##<<(character) column to group by. Must be a key in "groupKey" object
colorBy,	##<<(character) variable to colour by. Must be a key in "groupKey" object
oCol, 		##<<(character) group colour scheme
plotViewType, 	##<<(character) plot type for Data Track. See Gviz for options
plotType, 	##<<(character) currently hard-coded to "smoo2"
compEB, # T/F whether error bars should be computed
errb, # T/F whether error bars should be shown (must first be computed)
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
	if (colorBy=="(none)") { colgroups <- groupKey[["sampleName"]]; colorBy <- "sampleName" }
	else colgroups <- groupKey[[colorBy]]
} else {
	mygroups <- groupKey[[groupBy]]
	colgroups <- mygroups
}

# x/y limits, smoothing param
xlim <- c(min(outdat$start), max(outdat$start))  
if(whichYlim == "def"){
	tmp_mat <- as.matrix(outdat[,-c(1,2,grep("min",colnames(outdat)), grep("max",colnames(outdat)))])
	tmp <- na.omit(as.numeric(tmp_mat)); myYlim <- quantile(tmp,c(0.005,0.99)); 
} else myYlim <- customYlim
bw <- param2
xvals <- outdat$start

if (verbose) cat("\tset up xylim successfully\n")

# colours
colBase <- suppressWarnings(brewer.pal(length(colgroups)+1,oCol)[-1])
if (length(colBase) > length(colgroups)) colBase <- colBase[1:length(colgroups)]
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
	
	
	outdat <- outdat[,-(1:2)]
	
	# set colour and fill
	filldat <- colist
	if ("a_confint" %in% plotViewType){
		source("makeColorLighter.R")
		filldat <- sapply(colist, makeColorLighter)
	} 
	
	# smoothing if necessary
	vals <- outdat
	if (plotType=="smoo2") {
		vals <- smooMat(outdat,bw,xvals,1:ncol(outdat))
	}
	mcols(plot_GR) <- vals;
	
	# create data track
	cat("\t* Created data track\n")
	g <- groupKey[[groupBy]]; g <- g[which(g %in% mygroups)]
	cat(sprintf("Plot view type=%s\n", plotViewType))
	
	cat("going to get DataTrack\n")
	dTrack <- DataTrack(plot_GR,groups=factor(myfiles[,groupBy],levels=g),
		col=colist,fill=filldat, type=plotViewType,title="score")
	cat("Got dataTrack\n")
	
	# ######################################################
	# View mode 2: Individual samples
} else{ 
	# reorder based on group priority
  
	smpsGrps <- getSampleOrder(myfiles,groupKey)
	midx <- match(myfiles$sampleName, colnames(outdat))
	if (all.equal(colnames(outdat)[midx],myfiles$sampleName)!=T) { cat("sample names don't match output matrix column names. Check?")}
	outdat <- outdat[,midx]
	if (nrow(myfiles)==1) { outdat <- as.matrix(outdat); colnames(outdat) <- myfiles$sampleName }

	# assign sample colours
	colist <- colist[match(smpsGrps[,colorBy], names(colist))]
	
	if (plotType == "smoo2") outdat <- smooMat(outdat,param2, xvals,1:ncol(outdat))
	mcols(plot_GR) <- outdat
	dTrack <- DataTrack(plot_GR, legend=legd,groups=smpsGrps[,"sampleName"],
		aggregateGroups=F,type=plotViewType,
		col=colist)
} 
cat("setting display track\n")
displayPars(dTrack)  <- list(
	cex.title=1.5,background.title='cyan4', 
	baseline=0, col.baseline='gray20',lty.baseline=3,legend=legd,cex.axis=1.5,lwd=2,
	col.grid='gray50',lty.grid=3,cex=1,pch=20,
	grid=TRUE,ylim=myYlim
	)

eleft <- max(10,0.1 * diff(xlim))
finalTracks <- list(itrack,gtrack,dTrack)

sizes <- c(0.05,0.15,1.0)
cat("* About to get anno\n")
if (!is.null(selAnno)) {
	cat("\t* Collecting annotation tracks\n")
	annoTracks <- getAnnoTracks(configParams$annoConfig, selAnno,GRanges(setchrom,IRanges(xlim[1],xlim[2])),configParams$genomeName)
	finalTracks <- c(finalTracks, annoTracks$tracks)
	sizes <- c(sizes, annoTracks$sizes)
}

if (xlim[2]-xlim[1] > 3e5) gsize <- 0.05 else gsize <- 0.1
ttl <- sprintf("%s: %s, %1.1f-%1.1f Mb ; Groups: {%s}", plotTxt,setchrom, xlim[1]/1e6, xlim[2]/1e6, paste(mygroups,sep=",",collapse=","))
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
			if ((end(selGR)-start(selGR)) < 3e5) { 
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


