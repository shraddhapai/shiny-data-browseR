.. _install:
=======================
Install the Browse-R
=======================

This step only needs to be done once by the web administrator setting up the application. It does not have to be repeated by all endusers.
You will need access to a machine running a web server viewable to your users.


#. R dependencies:
	#.  Install shiny. When complete you should be able to run one of the built-in examples from the shiny package on your hosting machine.
	#.  Install Bioconductor and the following packages: 
		#. GenomicRanges
		#. Gviz
		#. Rsamtools
		#. any BioC annotation packages you may need. e.g. For gene annotation you may need TxDb.Mmusculus.UCSC.mm9.knownGene
#.  Get a copy of the Shiny BrowseR code from git: git clone ...
#. Move it to your shiny web directory. On most machines this is /var/shiny-server/www



