# ================================================================
#  Dockerfile — MR Pipeline CDC 1.0.0
#  Author : Nadeem Khan, Bioinformatician, INRS-CAFSB
# ================================================================

FROM rocker/r-ver:4.4.1

LABEL maintainer="Nadeem Khan <github.com/nkhan119>"
LABEL description="CDC 1.0.0 MR Pipeline — TwoSampleMR + MVMR + heterogeneity"
LABEL version="1.0.4"

# ── 1. System libraries ──────────────────────────────────────
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    pkg-config \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libgit2-dev \
    libssh2-1-dev \
    libxt-dev \
    libcairo2-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libzstd-dev \
    libgmp-dev \
    libmpfr-dev \
    python3 \
    python3-pip \
    python3-dev \
    wget \
    curl \
    unzip \
    pigz \
    tabix \
    fonts-liberation \
    ca-certificates \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# ── 2. plink2 + plink 1.9 ───────────────────────────────────
RUN wget -q https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_x86_64_20241206.zip -O /tmp/plink2.zip \
    && unzip -q /tmp/plink2.zip -d /usr/local/bin/ \
    && chmod +x /usr/local/bin/plink2 \
    && rm /tmp/plink2.zip \
    && wget -q https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20231211.zip -O /tmp/plink.zip \
    && unzip -q /tmp/plink.zip plink -d /usr/local/bin/ \
    && chmod +x /usr/local/bin/plink \
    && rm /tmp/plink.zip

# ── 3. Python packages ───────────────────────────────────────
RUN pip3 install --no-cache-dir \
    "numpy>=1.26" \
    "pandas>=2.1" \
    "matplotlib>=3.8" \
    "scipy>=1.11"

# ── 4. R CRAN layer — general utilities ──────────────────────
# Added 'parallelly' here to satisfy MVMR dependency requirements
RUN Rscript -e "\
  options(repos = c(CRAN = 'https://cloud.r-project.org')); \
  install.packages(c( \
    'remotes', 'BiocManager', 'data.table', 'R.utils', \
    'parallelly', 'dplyr', 'tidyr', 'Matrix', 'httr', \
    'httr2', 'jsonlite', 'curl', 'RCurl', 'xml2', \
    'openssl', 'Rcpp', 'RcppArmadillo', 'foreach', \
    'doParallel', 'MASS', 'lattice', 'survival', \
    'quantreg', 'mr.raps', 'cowplot', 'gridExtra', \
    'pbapply', 'psych', 'reshape2' \
  ), Ncpus = 4L) \
  "

# ── 5. R CRAN layer — gmp chain for MendelianRandomization ──
RUN Rscript -e "\
  options(repos = c(CRAN = 'https://cloud.r-project.org')); \
  chk <- function(p) { if (!requireNamespace(p, quietly = TRUE)) stop(paste('FAILED:', p)) }; \
  install.packages('gmp', Ncpus = 2L); chk('gmp'); \
  install.packages('arrangements', Ncpus = 2L); chk('arrangements'); \
  install.packages('iterpc', Ncpus = 2L); chk('iterpc'); \
  install.packages('MendelianRandomization', Ncpus = 2L); chk('MendelianRandomization'); \
  "

# ── 6. R Bioconductor ────────────────────────────────────────
RUN Rscript -e "\
  BiocManager::install(c('VariantAnnotation', 'GenomicRanges'), \
    ask = FALSE, update = FALSE) \
  "

# ── 7. R GitHub packages via codeload tarballs ───────────────
RUN Rscript -e "\
  options(repos = c(CRAN = 'https://cloud.r-project.org')); \
  \
  install_tarball <- function(url, pkg) { \
    message('Installing ', pkg, '...'); \
    tmp <- tempfile(fileext = '.tar.gz'); \
    download.file(url, tmp, quiet = FALSE, method = 'libcurl'); \
    # Using remotes::install_local is more robust for GitHub tarballs
    remotes::install_local(tmp, dependencies = TRUE, upgrade = 'never', INSTALL_opts = '--no-multiarch'); \
    unlink(tmp); \
    if (!requireNamespace(pkg, quietly = TRUE)) stop(paste('FAILED:', pkg)); \
  }; \
  \
  install_tarball('https://codeload.github.com/mrcieu/ieugwasr/tar.gz/refs/heads/master', 'ieugwasr'); \
  install_tarball('https://codeload.github.com/rondolab/MR-PRESSO/tar.gz/refs/heads/master', 'MRPRESSO'); \
  install_tarball('https://codeload.github.com/WSpiller/RadialMR/tar.gz/refs/heads/master', 'RadialMR'); \
  install_tarball('https://codeload.github.com/gqi/MRMix/tar.gz/refs/heads/master', 'MRMix'); \
  install_tarball('https://codeload.github.com/MRCIEU/TwoSampleMR/tar.gz/refs/heads/master', 'TwoSampleMR'); \
  install_tarball('https://codeload.github.com/WSpiller/MVMR/tar.gz/refs/heads/master', 'MVMR'); \
  "

# ── 8. Final verification ────────────────────────────────────
RUN Rscript -e "\
  pkgs <- c('TwoSampleMR', 'MendelianRandomization', 'MVMR', 'ieugwasr', 'MRPRESSO', 'RadialMR', 'MRMix'); \
  for (p in pkgs) { if (!requireNamespace(p, quietly = TRUE)) stop(paste('Verification failed for:', p)) }; \
  cat('All R packages verified successfully.\n'); \
  "

WORKDIR /pipeline
CMD ["/bin/bash"]
