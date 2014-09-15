# BSseq.R
# fetches data for BSseq objects

fetchData_base <- function
( 
session,
pheno, 		##<<(data.frame) phenotype matrix
selRange, ##<<(GRanges) range being viewed on browser
bin_GR,  ##<<(GRanges) ranges of individual data bins
numBins,	##<<(integer) num. bins
aggFUN=mean,	##<<(function) aggregating function
myParams  ##<<(list) format-specific params. For bigwig this is empty and does nothing.
# CHROM_POS - columm # of sequence name in tabix file
# START_POS - columm # of position start in tabix file
# END_POS - column # of position end in tabix file
# STRAND_POS - column # of strand
# M_POS - column # of num. methylated cytosine (M)
# COV_POS - column # of position coverage
# minCov - minimum coverage to use
) {

  verbose <- FALSE
	# initialize progress bar
	updateProgressBar(session, inputId="load_pBar", visible=TRUE,value=0,color="warning",striped=TRUE,animate=F)

  
myVals <- list(CHROM_POS=1,START_POS=2,END_POS=2,STRAND_POS=3,M_POS=4, COV_POS=5,minCov=4)
for (nm in names(myParams)) {
   myVals[[nm]] <- as.integer(myParams[[nm]])
}

tabixFiles <- pheno$bigDataURL
tLen <- length(tabixFiles)
out_score <- matrix(nrow=length(bin_GR), ncol=length(tabixFiles))

ctr <- 1
for(t in tabixFiles) { 
	# read tabix.
	if (verbose) cat("\t scanning\n")
	t0 <- system.time(tbx <- scanTabix(TabixFile(t),param=bin_GR))
	if (verbose) print(t0)
	cat(sprintf("Got %i records\n", length(tbx)))

	if (verbose) cat("\t computing %M\n")
	M_POS <- myVals[["M_POS"]]; COV_POS <- myVals[["COV_POS"]]
	sc <- mclapply(1:length(tbx), function(curr_idx) { # process ranges in parallel
#	sc <- sapply(1:length(tbx),function(curr_idx) { # serial loop -- for debugging
		cur <- tbx[[curr_idx]]
		if (length(cur) < 1) return(NA)

		# convert records to data.frame.
		nc <- length(strsplit(cur[1],"\t")[[1]])
		rec <- strsplit(cur,"\t")
		rec <- as.data.frame(matrix(unlist(rec),byrow=T, ncol=nc),stringsAsFactors=F)
		rec[,M_POS]<- as.integer(rec[,M_POS])
		rec[,COV_POS] <- as.integer(rec[,COV_POS])

		idx <- which(rec[,COV_POS] >= myVals[["minCov"]])
		if (!any(idx)) return(NA)

		# compute score
		COV <- sum(rec[idx,COV_POS]); M <- sum(rec[idx,M_POS]);
		pctM <- (M/COV)*100
		#if (verbose) cat(sprintf("M = %i ;  Cov=%i ; %%M = %1.1f\n", M, COV, pctM))
		return(pctM)
	},mc.cores=4)
#	}) # end sapply - for debugging
	tryCatch({
		#out_score[,ctr:(ctr+(outcol-1))] <- matrix(unlist(sc),byrow=T,ncol=outcol); 
    out_score[,ctr] <- unlist(sc)
    ctr <- ctr+1

	updateProgressBar(session, inputId="load_pBar", value=round((ctr/tLen)*100),color="warning",
							striped=TRUE,animate=FALSE)
	
    #ctr <- ctr+outcol
	}, error=function(ex) {
		print(ex); browser()
	})
	} 
	updateProgressBar(session, inputId="load_pBar", value=100,color="success",
							striped=FALSE,animate=FALSE)
	

colnames(out_score) <- pheno$sampleName
return(out_score)
}

getSeqinfo <- function(
  f ##<<(character) a representative file of the dataset. All files are assumed to have the same sequences represented
) {
  hdr <- headerTabix(f);
  return(hdr$seqnames)
}
