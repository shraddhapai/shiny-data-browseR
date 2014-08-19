#Scripts related to data processing

listConfig <- function# list all config files available
(){
cat(sprintf("Currently I'm in %s",getwd()))
  indir <- scan("config_location.txt",what="character")
cat(sprintf("Root directory is: %s\n", indir))
fList <- dir(path=indir,pattern="_config.txt")
configSet <- list()
for (f in fList) { 
  cat(sprintf("Config file: %s\n", f))
	tmp <- getSettings(sprintf("%s/%s", indir,f)); 
	print(tmp$name)
	configSet[[tmp$name]] <- tmp 
}

return(configSet)
### (list of lists) key: name of dataset, value: config key-value pairs
}

getSettings <- function# parse config file
(configFile="config.txt", printMe=F
) {
  dat <- read.delim(configFile,sep="\t",header=F,as.is=T,comment.char="#")
  configParams <- list()
  for (k in 1:nrow(dat)) {
    configParams[[dat[k,1]]] <- dat[k,2]
    if (printMe) cat(sprintf("%s\t%s\n", dat[k,1],dat[k,2]))
	  if (dat[k,1]=="groupCols") configParams[[dat[k,1]]] <- unlist(strsplit(dat[k,2],","))
    }
if (!file.exists(configParams$datasetConfig)) stop("Invalid dataset-specific config file. Add 'datasetConfig' key to config.txt. \n")
if (!file.exists(configParams$chromSizes)) stop("chromSizes file not specified\n");
if (!file.exists(configParams$groupOrder)) stop("Invalid group order file. Add valid 'groupOrder' key to config.txt\n")
if (!"groupCols" %in% names(configParams)) stop("Specify groupCols parameter so I know what columns contain sample groups\n")
if (!"datatype" %in% names(configParams)) stop("Specify 'datatype' parameter so you datatype can be correctly handled.")

# move datatype-specific vars to different list
data_idx <- grep(sprintf("^%s__", configParams$datatype), names(configParams))
datatypeParams <- list()
if (any(data_idx)) {
  for (nm in names(configParams)[data_idx]) {
      nm_new <- sub(sprintf("^%s__",configParams$datatype),"", nm)
      datatypeParams[[nm_new]] <- configParams[[nm]]
      configParams[[nm]] <- NULL
  }
}
                 
pheno <- read.delim(configParams$datasetConfig, header=T, as.is=T, sep="\t")
groupInfo <- read.delim(configParams$groupOrder,sep="\t",header=T,as.is=T,comment.char="#")
chromsize <- read.delim(configParams$chromSizes, header=F, as.is=T, sep="\t"); colnames(chromsize) <- c("chrom", "size")

groupKey <- list()
for (k in 1:nrow(groupInfo)) {
	groupKey[[groupInfo[k,1]]] <- unlist(strsplit(groupInfo[k,2],","))
}
for (k in setdiff(configParams$groupCols, names(groupKey))) { # for groups without specified order
	tmp <- sort(unique(pheno[,k])); ; # just put them in alphanumeric order
	groupKey[[k]] <- tmp
}

 return(list(name=configParams$name, allDat=pheno, groupKey=groupKey, 
             chromsize=chromsize, configParams=configParams,
             datatypeParams=datatypeParams))
}

fetchData <- function# Extracts data for selected range from all BW files, returns as a data.frame
(
session, 
myfiles,	##<<(data.frame) sample phenotype df
selRange,	##<<(GRanges) genomic interval to extract
numBins,	##<<(integer) num. bins
datatype,	##<<(character) 'datatype' variable from dataset config file
verbose=F,
FUN=mean, 
datatypeParams ##<<(list) params specific to this datatype
) {
	# prepare output matrix. coords + score from first file
	verbose <- TRUE
	if (verbose) cat("In fetchData\n")
	x1 <- start(selRange)[1]; x2 <- end(selRange)[1]
        binSize <- floor((x2-x1)/numBins)
	binStarts <- seq(x1,x2-binSize, binSize); binEnds <- c(binStarts[-1]-1, x2)
	bin_GR <- GRanges(seqnames(selRange)[1], IRanges(binStarts,binEnds))

	parserFile <- sprintf("data_types/%s.R", datatype)
		if (!file.exists(parserFile)) {
			stop(sprintf("Parser file for datatype %s not found!", datatype))
		}

		source(parserFile);
		tryCatch({
			out <- fetchData_base(session=session, pheno=myfiles, selRange=selRange, bin_GR=bin_GR, 
                            numBins=numBins, aggFUN=mean, myParams=datatypeParams)
		},error=function() {
			cat(sprintf("Error in fetchData() function for datatype: %s", datatype))
			print(ex)
		}, finally={} 
		)
		bed <- data.frame(chrom=seqnames(selRange)[1], start=binStarts, end=binEnds)
		if (verbose) cat("fetchData::about to return\n")
return(list(coords=bed,values=out))
}

computeAverages <- function# computes group mean and optionally variation
(
myfiles,	##<<(data.frame) sample pheno table
mygroups,	##<<(character) vector of group names
bed,		##<<(data.frame) coordinates
alldat,		##<<(matrix) sample-wise scores
compEB,		##<<(logical) if T compute CI; else don't
groupBy,		##<<(character) which grouping to use?
verbose=F
) {
t0 <- system.time(
grp.summ <- foreach(i = 1:length(mygroups), .combine=cbind, .inorder=F) %dopar% {
          g <- mygroups[i]
          thisgroup <- myfiles$sampleName[which(myfiles[,groupBy] == g)]
          idx <- which(colnames(alldat) %in% thisgroup)
		  if (verbose) cat(sprintf("\t %s: %i samples\n",g, length(idx))) 
          
          ## will get error messages since some rows do not have any non-NA values (at the group level)
          ## this will produce stats of NA, and will not affect future analysis
          ## don't want to eliminate here in order to keep all data together
          
          if(length(idx)==1){
            grp.mean <- alldat[,idx]
            grp.max <- alldat[,idx]
            grp.min <- alldat[,idx]
          }else{
            grp.mean <- suppressWarnings(rowMeans(alldat[,idx], na.rm=T))

            if(compEB == TRUE){
				stop("computing EB is currently disabled")
            }else{ # to avoid error messages since expecting something
              grp.min <- grp.mean 
              grp.max <- grp.mean
            }
          }
          tmp <- data.frame(grp.mean, grp.max, grp.min)
          colnames(tmp) <- paste(g, c("mean", "max", "min"), sep=".")
          
          return(tmp)
        })
	if (verbose) { cat("Grouping time:\n"); print(t0) }

return(grp.summ)
}

baselineSamps <- function
(
alldat,			##<<(data.frame) matrix of sample-wise intensities
grp.summ,		##<<(data.frame) group mean,min,max for all groups
mygroups,		##<<(character) group names
whichMetric, 	##<<(character) which metric to compute for sample-wise data;  options are [normal|logratio_all|logratio_baseline]
whichBaseline,	##<<(character) group for baseline (for logratio_baseline option)
logMe			##<<(logical) log-transform?
) {
if (whichMetric=="ratio_baseline") 
	mybase <- grp.summ[,sprintf("%s.mean",whichBaseline)]
else 
	mybase <- rowMeans(alldat,na.rm=T)

alldat <- alldat/mybase
if (logMe) alldat <- log2(1+alldat);


for (g in mygroups) {
	x <- sprintf("%s.mean",g)
	grp.summ[,x] <- grp.summ[,x]/mybase
	if (logMe) grp.summ[,x] <- log2(1+grp.summ[x])
}

return(list(alldat=alldat,grp.summ=grp.summ))
}

smooMat <- function
(
dat,
bw,
xpoints,
col2run
) {
 smooOut <- foreach(i=col2run, .combine=cbind, .inorder=T) %dopar% {
 	mycol <- dat[,i]
	rm_idx <- which(is.na(mycol))
    if (!any(rm_idx)) ksmoo <- ksmooth(x=xpoints, y=mycol,bandwidth=bw,x.points=xpoints)
    else ksmoo <- ksmooth(x=xpoints[-rm_idx], y=mycol[-rm_idx],bandwidth=bw, x.points=xpoints)
    return(ksmoo[["y"]])
	}
if (is.null(ncol(smooOut))) smooOut <- as.matrix(smooOut);
colnames(smooOut) <- colnames(dat)[col2run]
return(smooOut)
}

getSampleOrder <- function# assign sample priority based on provided group ordering
(smpsGrps, ##<<(data.frame) sample group info. each row is sample. 
groupKey
) {
myscore <- rep(0,nrow(smpsGrps))
gnames <- names(groupKey)
for (k in 1:nrow(smpsGrps)) {
	for (g in 1:length(groupKey)) {
    cat(sprintf("k=%s\ng=%s\n\n", paste(smpsGrps[k,],collapse=",",sep=","),gnames[g]))
    idx <- which(groupKey[[g]]==smpsGrps[k,gnames[g]])
    if (!any(idx)) {
      
     s1 <- sprintf("****************\nERROR: Group member %s is not a listed member of group %s\n", smpsGrps[k,gnames[g]], gnames[g])
     s1 <- sprintf("%s Offending sample:\n%s\n", s1, print(smpsGrps[k,],collapse=",",sep=","))
     s1 <- sprintf("%s Suggestions: update group_order.txt to include this entry or change it to a valid member.", s1) 
     stop(s1)
    }
		myscore[k] <- myscore[k] + (10^(g-1))*which(groupKey[[g]] == smpsGrps[k,gnames[g]])
	}
}
smpsGrps <- smpsGrps[order(myscore),]
}


