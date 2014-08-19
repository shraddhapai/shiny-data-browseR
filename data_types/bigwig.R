# bigwig.R
# fetch data for bigwig data format

fetchData_base <- function
(
session, 
pheno, 		##<<(data.frame) phenotype matrix
selRange, 	##<<(GRanges) range being viewed on browser [start,end] - length 1
bin_GR,		##<<(GRanges) ranges of individual data bins
numBins,	##<<(integer) num. bins
aggFUN=mean,	##<<(function) aggregating function
myParams  ##<<(list) format-specific params. For bigwig this is empty and does nothing.
) {
	verbose <-TRUE 
	# initialize progress bar
	updateProgressBar(session, inputId="load_pBar", visible=TRUE,value=0,color="warning",striped=TRUE,animate=F)

	 # now for all other files:
	x1 <- start(selRange)[1]; x2 <- end(selRange)[1]
        binSize <- floor((x2-x1)/numBins)
	cat("about to get samples\n")
	alldat <- foreach(s=1:nrow(pheno), .combine=cbind, .inorder=T) %dopar% {
          gr <- import.bw(pheno$bigDataURL[s],asRangedData=F,which=selRange)
		  if (verbose) cat(sprintf("%s: %i records\n", s, length(gr)))
		  
		  # bin data
		  binNum <- floor((start(gr)-x1)/binSize)+1
		  too_far <- which(binNum > length(bin_GR)); if (any(too_far)) binNum[too_far] <- numBins
		  x <- aggregate(gr$score, by=list(binNum=binNum), FUN=aggFUN); 

		  out <- numeric(length=length(bin_GR))+NA; out[x$binNum] <- x$x

		  updateProgressBar(session, inputId="load_pBar", value=round((s/nrow(pheno))*100),color="warning",
							striped=TRUE,animate=FALSE)

          return(out)
     }
	if (nrow(pheno)==1) alldat <- as.matrix(alldat)
	colnames(alldat)<- pheno$sampleName 
	updateProgressBar(session, inputId="load_pBar", value=100,color="success",
							striped=FALSE,animate=FALSE)


return(alldat)
### (matrix) sample-wise values. Row order should correspond to coords in bin_GR and column order to samples.
}

getSeqinfo <- function(
  f ##<<(character) a representative file of the dataset. All files are assumed to have the same sequences represented
  ) {
  hdr <- seqinfo(BigWigFile(f))
  return(hdr@seqnames)
}
