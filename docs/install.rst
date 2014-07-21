.. _install:

.. toctree::
	:maxdepth: 2
	:numbered:

=======================
Install the Browse-R
=======================

This step only needs to be done once by the web administrator setting up the application. It does not have to be repeated by all endusers.
You will need access to a machine running a web server viewable to your users.


#. Install R (>=3.0.1)
#. Install the following R packages. Some are from the general CRAN repository, while others are from Bioconductor.
	#. shiny (CRAN); (>= 0.9.1)
	#. RColorBrewer (CRAN)
	#. Biobase (BioC; >=2.22.0)
	#. Gviz (BioC; >=1.6.1)
	#. GenomicRanges (BioC; >=1.14)
	#. Rsamtools (BioC; >=1.14)
	#. (optional) any BioC annotation packages you may need. e.g. For gene annotation you may need TxDb.Mmusculus.UCSC.mm9.knownGene
#.  Get a copy of the Shiny BrowseR code from git: git clone ...
#. Move it to your shiny web directory. On most machines this is /var/shiny-server/www

========================
Add data sources
========================
To add a new dataset to the browser, you need to setup the following:

#. *dataset config file*: A dataset-specific config file containing all basic informatino about the dataset
#. *phenotype matrix*: A tab-delimited text file containing a table of sample-wise metadata and location of sample-wise data file. See section X.XX for details.
#. *data*: The data itself - e.g. sample-wise epigenetic measures. Currently this data is expected to be in bigwig format; future releases of the browser are expected to support other file formats.
#. *group order*: A file indicating the order in which elements of a group are assigned colours in the browser. For example, a timeseries experiment may find it useful to cause the first timepoint to be plotted as first in the series, the second timepoint as second in a series. etc., rather than rely on alphanumeric ordering of the filenames. This is the file where the relative ordering of category members is specified.

**Note: The browser expects all dataset config files to be located in a single directory.** 
The location of this directory is specified in "config_location.txt", which is located in the directory from which the app is run: e.g. */var/shiny-server/www/Shiny_BrowseR*

Example directory structure for datasets
----------------------------------------------------------

.. _add-data-exampledir:

*TODO Mention no constraint to have data in this structure, as links are provided for data files in phenotype matrix and config files can be symlinked.*
In this example we have two datasets, one showing epigenetic dynamics in mouse brain development and the other, for human. Here, :code:`config_location.txt` points to a directory which also contains all browser-related metadata and the data files themselves::

	<config_location_points_here>
		 |
		 |---- mouseBrainDev_config.txt *symlink*
		 |---- humanBrainDNAMethylome_config.txt *symlink*
		 |---- mouseBrainDev/
		 		|---- config.txt
				|---- pheno.txt
				|---- group_order.txt
				|---- data
					  |---- GSM123456_mm9.bw
					  |---- GSM654321_mm9.bw
					  | ...
					  |---- GSM9999_mm9.bw
		 |---- humanBrainDNAMethylome/
		 		|---- config.txt
				|---- pheno.txt
				|---- group_order.txt
				|---- data
					  |----
					  |---- PCW4.bw
					  |---- midGestation_Cortex.bw
					  | ...
					  |---- PostPuberty_AnteriorCingulateCortex.bw
	     |----anno


*TODO: Turn above into a graphic?*

Dataset config file
-----------------------

.. _add-data-config:

This file contains metadata about the dataset in general. It is a tab-delimited file with two columns: a key (controlled word) and value.
Currently, all these fields are required:
* **TODO: Fill this in**

:ref:`Back to top <install>`

Phenotype matrix
------------------------

.. _add-data-pheno:

This tab-delimited file contains sample-wise metadata, including locations of data files. Each row should contain data for one sample, and each column should contain a unique type of metadata. Column order is unimportant to the browser.
The browser expects the following columns, named exactly in this way:

* sampleName *TODO or is it SampleID? Look it up!*
* bigDataURL: Location of data file
* all grouping columns as described in the :ref:`grouping order <add-data-grouping>` file.

Groups and grouping order
---------------------------

.. _add-data-grouping:

The :code:`group_order.txt` file is a tab-delimited file containing a table of two columns:
#. groupID: Group name, must match a column name in the phenotype matrix
#. groupOrder: Order in which group members must be shown. Comma-separated collection of values. **All values for a given group must be specified here. The browser will return an error if any additional group members are found in the phenotype table but are not listed here.

In addition to these groups, the browser allows a non-grouping option - i.e. viewing sample-specific data - with "Grouping: (none)". 

As an example::

	groupID	groupOrder
	Tissue	Brain,Sperm
	Diagnosis	Control,Schizophrenia,Bipolar disorder
	TimeOfSampling	Before_Treatment,During_Treatment,After_Treatment

For this dataset, the browser would show 4 grouping options: Tissue, Diagnosis, TimeOfSampling, (none).

 :ref:`Back to top <install>`


========================
Add annotation sources
========================

This is the directory structure for annotation sources::

		 	   |------ hg19
			   		   |------ cpgIslandExt.txt
					   |------ cytoBandIdeo.txt
					   |------ TxDb.Hsapiens.UCSC.hg19.refGene.sqlite
			   |------ mm9
			   		   |------ cpgIslandExt.txt
					   |------ LAD_NPC_mm9.txt

This is normal text again.

:ref:`Back to top <install>`
