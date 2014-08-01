.._add_datatype: 

============================
Extending Browse-R capabilities
============================

Add a datatype handler
-----------------------

The Shiny Browse-R is easily customizable to show data from a new platform; all that is required is the addition of a platform-specific .R file in the code directory. For illustration, let us suppose we want to add a new datatype of 450K microarrays.

These steps will allow a new 'datatype' value to be specified in a :ref:`dataset config <add-data-config>` file. In our example, let that datatype be "FourFiftyK".

#. Create FourFiftyK.R in the ``data_types`` directory; the latter is at the same level with ui.R and server.R)
#. In FourFiftyK.R, create an R function, fetchData_base() with the following signature::

	fetchData_base <- function
	(
	pheno, 		##<<(data.frame) phenotype matrix
	selRange, 	##<<(GRanges) range being viewed on browser [start,end] - length 1
	bin_GR,		##<<(GRanges) ranges of individual data bins
	numBins,	##<<(integer) num. bins
	aggFUN=mean	##<<(function) aggregating function
	) {
		# fetch code goes here.
	
	### (list with two keys):
	### 1) coords: data.frame with three columns corresonding to the chromosome, start and end coordinate
	### 2) values: sample-wise values. Row order should correspond to coords and column order to samples.
}

Current datatypes
-------------------
bigwig
^^^^^^
Used for any datatype which can be represented with a single column containing a continuous value.

BSseq
^^^^^
Used for bisulfite-seq data. Usually has two columns, M and COV, which are combined into a %methylation over an arbitrary genomic interval (e.g. in 2Kb bins, or over a gene).

Example view of input file: 

Custom columns are: 
# CHROM_POS - columm # of sequence name in tabix file
# START_POS - columm # of position start in tabix file
# END_POS - column # of position end in tabix file
# STRAND_POS - column # of strand
# M_POS - column # of num. methylated cytosine (M)
# COV_POS - column # of position coverage
# minCov - minimum coverage to use

Computation is: (M/COV)

Future version will incorporate non-conversion rate subtraction capability.
