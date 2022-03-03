FROM rocker/verse:4.0.1

# system dependencies
RUN export DEBIAN_FRONTEND=noninteractive; apt-get -y update \
  && apt-get install -y gdal-bin=3.0.4+dfsg-1build3 \
	libgdal-dev=3.0.4+dfsg-1build3 \
	libgeos-dev=3.8.0-1build1 \
	libgeos++-dev=3.8.0-1build1 \
	libudunits2-dev=2.2.26-5 \
	make=4.2.1-1.2 \
	pandoc=2.5-3build2 \
	pandoc-citeproc=0.15.0.1-1build4

# working directory name
WORKDIR /home/rstudio/hkl-mirs

# install renv R package
ENV RENV_VERSION 0.11.0-6
RUN R -e "install.packages('emotes', repos = c(CRAN = 'https://cran.r-project.org'))"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

# Restore renv cache
COPY renv.lock renv.lock
RUN chown -R rstudio . \
 && sudo -u rstudio Rscript -e 'renv::consent(provided = TRUE)' \
 && sudo -u rstudio Rscript -e 'Sys.setenv(R_INSTALL_STAGED = FALSE)' \
 && sudo -u rstudio R -e 'renv::restore(repos = c(RSPM = "https://packagemanager.rstudio.com/all/latest", CRAN = "https://cran.r-project.org/"))'

# install compendium
RUN sudo -u rstudio Rscript -e 'remotes::install_github("henningte/hklmirs")'

# labels
LABEL maintainer = "Henning Teickner <henning.teickner@uni-muenster.de>" \
  org.opencontainers.image.authors = "Henning Teickner <henning.teickner@uni-muenster.de>" \
  author.orcid = "0000-0002-3993-1182" \
  org.opencontainers.image.version = "0.0.0.9000" \
  org.opencontainers.image.licenses = "GPL-3"

# instructions

# to build the image, navigate to the directory with the Dockerfile and run (the images has a size of ~7 Gb):
# docker build -t hklmirs:0.0.0.9000 .

# to run the image in a container, do:
# docker run --rm -e PASSWORD=hkl -p 8787:8787 -v $(pwd):/home/rstudio --user rstudio hklmirs:0.0.0.9000

# to reproduce the analyses, run (takes about 12 h and occupies additional 2 Gb disk space)
# docker run --rm -v $(pwd):/home/rstudio hklmirs:0.0.0.9000 R -e 'setwd("/home/rstudio/hkl-mirs"); rmarkdown::render("analysis/paper/001-paper-main.Rmd");'


