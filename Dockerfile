# Base image: R 4.4.2 (Based on Ubuntu 24.04 LTS)
FROM rocker/r-ver:4.4.2

# Install Ubuntu System Dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libglpk-dev \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    libudunits2-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libhdf5-dev \
    libgit2-dev \
    cmake \
    gfortran \
    patch \
    && rm -rf /var/lib/apt/lists/*

# Configure P3M
RUN echo "options(repos = c(CRAN = 'https://packagemanager.posit.co/cran/__linux__/noble/latest'))" >> /usr/local/lib/R/etc/Rprofile.site

ENV R_LIBS_SITE=/usr/local/lib/R/site-library

RUN R -e "install.packages(c('remotes','BiocManager'))"
RUN R -e "BiocManager::install(version = '3.20', ask = FALSE)"

# Install CRAN Packages
# Combined list of all your CRAN requirements
RUN R -e "remotes::install_version('gtools', version = '3.9.5')" \
    && R -e "remotes::install_version('viridis', version = '0.6.5')" \
    && R -e "remotes::install_version('colorspace', version = '2.1-1')" \
    && R -e "remotes::install_version('ggrepel', version = '0.9.6')" \
    && R -e "remotes::install_version('tidyverse', version = '2.0.0')" \
    && R -e "remotes::install_version('matrixStats', version = '1.5.0')" \
    && R -e "remotes::install_version('SeuratObject', version = '5.0.2')" \
    && R -e "remotes::install_version('Seurat', version = '5.2.1')" \
    && R -e "install.packages(c('quantmod', 'pracma', 'hdf5r', 'zoo', 'scales', 'ggsignif', 'ggtext', 'devtools', 'Matrix', 'data.table'))"

RUN R -e 'install.packages( \
  "https://cran.r-project.org/src/contrib/Archive/grr/grr_0.9.5.tar.gz", \
  repos = NULL, \
  type = "source" \
)'

# Install Bioconductor Dependencies
RUN R -e "BiocManager::install(c( \
    'SingleCellExperiment', \
    'SummarizedExperiment', \
    'GenomicRanges', \
    'GenomeInfoDb', \
    'IRanges', \
    'S4Vectors', \
    'MatrixGenerics', \
    'Biobase', \
    'BiocGenerics', \
    'limma', \
    'DelayedArray', \
    'DelayedMatrixStats', \
    'batchelor', \
    'biomaRt', \
    'rhdf5', 'HDF5Array' \
    ), update = FALSE, ask = FALSE)"

# 6. Install Monocle3
RUN R -e "remotes::install_github('cole-trapnell-lab/monocle3', upgrade='never')"
#RUN R -e "library(monocle3)"

# 7. Setup Working Directory
WORKDIR /analysis
COPY . /analysis

# 8. Entrypoint
ENTRYPOINT []
CMD ["/bin/bash"]