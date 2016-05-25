FROM ubuntu:14.04

MAINTAINER Etienne Thevenot (etienne.thevenot@cea.fr)

# Setup package repos
RUN echo "deb http://mirrors.ebi.ac.uk/CRAN/bin/linux/ubuntu trusty/" >> /etc/apt/sources.list
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9

# Update and upgrade system
RUN apt-get update
RUN apt-get -y upgrade

# Install R and other needed packages
RUN apt-get -y install r-base
RUN R -e "install.packages('ropls', lib='/usr/lib/R/library', dependencies = TRUE, repos='http://mirrors.ebi.ac.uk/CRAN')"

# Clone tool
RUN apt-get -y install git
RUN git clone -b docker https://github.com/workflow4metabolomics/multivariate /files/multivariate

# Clean up
RUN apt-get clean && apt-get autoremove -y && rm -rf /var/lib/{apt,dpkg,cache,log}/ /tmp/* /var/tmp/*

# Define Entry point script
ENTRYPOINT ["/files/multivariate/multivariate_wrapper.R"]
