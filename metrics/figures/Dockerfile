FROM rocker/r-ver:4

RUN apt-get update && apt-get install -y libcurl4-openssl-dev libxml2-dev libssl-dev 
ENV SRC /usr/local/src
ENV BIN /usr/local/bin

WORKDIR $SRC


#RUN Rscript -e 'Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")'
RUN Rscript -e 'install.packages("remotes",dependences=TRUE); remotes::install_github("katiejolly/nationalparkcolors")'
RUN Rscript -e 'install.packages("readxl",dependencies=TRUE)'
RUN Rscript -e 'install.packages("tidyr",dependencies=TRUE)'
RUN Rscript -e 'install.packages("plyr",dependencies=TRUE)'
RUN Rscript -e 'install.packages("dplyr",dependencies=TRUE)'
RUN Rscript -e 'install.packages("argparser",dependencies=TRUE)'
RUN Rscript -e 'install.packages("cowplot",dependencies=TRUE)'
RUN Rscript -e 'install.packages("viridis",dependencies=TRUE)'
RUN Rscript -e 'install.packages("reshape2",dependencies=TRUE)'
RUN Rscript -e 'install.packages("htmlwidgets",dependencies=TRUE)'
RUN Rscript -e 'install.packages("webshot",dependencies=TRUE); webshot::install_phantomjs()'
RUN Rscript -e 'install.packages("plotly",dependencies=TRUE)'


COPY combine_results.R /bin
COPY compare_sigs.R /bin

RUN chmod 777 /bin/combine_results.R
RUN chmod 777 /bin/compare_sigs.R

VOLUME ['/tmp']
