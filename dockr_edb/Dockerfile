FROM wch1/r-shiny-server

# Thanks to Winston Chang of Rstudio,Inc. for providing the 
# dockerfile for r-shiny-server
# https://github.com/wch/dockr/blob/master/r-shiny-server/Dockerfile 

MAINTAINER Shraddha Pai "Shraddha.Pai@camh.ca"

# so we can see debug messages in R
ENV LC_ALL C

# Needed for BioC rtracklayer package
RUN apt-get install -y \
    libxml2-dev \
    git

# Get EDB repo + devel version of Gviz + EDB demo data
ADD EDB_demodata.tar.gz /EDB_demodata.tar.gz   
# TODO: move to wget when issue is fixed
# wget http://epidatabrowser.camh.ca:3838/demo_datasets/EDB_demodata.tar.gz && \

RUN git clone https://github.com/shraddhapai/shiny-data-browseR.git 
    tar xvfz EDB_demodata.tar.gz && \
    mv EDB_demodata /home/docker/. && \
    wget http://www.bioconductor.org/packages/devel/bioc/src/contrib/Gviz_1.11.2.tar.gz 

# install R + BioConductor dependencies
RUN R -e "install.packages(c('shinyBS','doMC','RColorBrewer', 'matrixStats'), repos='http://cran.rstudio.com/'); " && \
    R -e "source('http://bioconductor.org/biocLite.R'); biocLite(c('Biobase','GenomicRanges','rtracklayer','GenomicFeatures','biovizBase','Rsamtools'));" && \
    R -e "install.packages('Gviz_1.11.2.tar.gz')"

# move EDB code to server directory and create file pointing to demo data location
RUN mkdir /srv/shiny-server/EDB && \
    cp /shiny-data-browseR/*.R /srv/shiny-server/EDB/. && \
    cp -r /shiny-data-browseR/www /srv/shiny-server/EDB/. && \
    cp -r /shiny-data-browseR/data_types /srv/shiny-server/EDB/. && \
    echo "/home/docker/EDB_demodata/" > /srv/shiny-server/EDB/config_location.txt

# start the server
RUN touch /srv/shiny-server/restart.txt
