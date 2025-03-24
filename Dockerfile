# Use the official Conda base image
FROM continuumio/miniconda3:4.12.0

# Set the working directory
WORKDIR /app

# Install system dependencies (Git and R-related dependencies)
RUN apt-get update && apt-get install -y \
    git \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    build-essential \ 
    && rm -rf /var/lib/apt/lists/*
    
# update all packages to conda-forge channel    
RUN conda upgrade -c conda-forge --all -y

# Install R and related packages
RUN conda install -c conda-forge -c bioconda \
    r-base=4.4.2 \
    r-r.utils \
    r-vcfr \
    r-devtools \
    bioconductor-rsamtools \
    bioconductor-genomicfeatures \
    bioconductor-deseq2 \
    bioconductor-qvalue \
    bioconductor-bsgenome.hsapiens.ucsc.hg38 \
    bioconductor-txdbmaker \
    bioconductor-guitar \
    bioconductor-chippeakanno \
    bioconductor-enhancedvolcano \
    bioconductor-clusterprofiler \
    bioconductor-org.hs.eg.db \
    -y

# # Install R packages from GitHub 
RUN R -e "devtools::install_github(c('scottzijiezhang/MeRIPtools', 'scottzijiezhang/RADAR', 'ZW-xjtlu/exomePeak2'), build_vignettes = FALSE, upgrade_dependencies = FALSE)"

# Set the default command 
# CMD ["Rscript", "--vanilla", "/app/script.R"]
CMD ["R"]