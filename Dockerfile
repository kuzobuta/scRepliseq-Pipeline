FROM ubuntu:16.04
LABEL maintainer "Hisashi Miura <hisashi.miura[at]riken.jp>"

RUN apt update -y && DEBIAN_FRONTEND=noninteractive apt install -y \
    apt-utils \
    bzip2 \
    gcc \
    git \
    less \
    dialog \
    libncurses-dev \
    make \
    time \
    unzip \
    vim \
    wget \
    zlib1g-dev \
    software-properties-common \
    apt-transport-https \
    language-pack-en \
    libxml2-dev \
    autoconf

RUN apt update -y && apt install -y \
    libcurl4-gnutls-dev \
    libssl-dev

RUN apt update -y && apt install -y bedtools=2.25.0-1

##Install R version 3.4.4 (2018-03-15) -- "Someone to Lean On"
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu xenial/'
RUN apt update -y && apt install -y r-base=3.4.4-1xenial0

##Install "AneuFinder_1.2.1"
RUN R -e 'install.packages("BiocInstaller", repos="http://bioconductor.org/packages/3.4/bioc")'
RUN R -e 'install.packages("ggplot2"); install.packages("gtable"); install.packages("plyr")'
RUN R -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/cowplot/cowplot_0.9.4.tar.gz", repos=NULL, method="libcurl", dependencies=TRUE)'
RUN R -e 'install.packages("bitops")'
RUN R -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/caTools/caTools_1.17.1.2.tar.gz", type="source", repos=NULL, dependencies=TRUE)'
RUN R -e 'library(BiocInstaller); biocLite("AneuFinder")'
RUN R -e 'install.packages("pracma"); install.packages("optparse"); install.packages("XML"); install.packages("zoo")'

WORKDIR /usr/local/bin
COPY Download_Programs.sh .
RUN . Download_Programs.sh

## installing Miniconda version 4.6.14
RUN wget https://repo.continuum.io/miniconda/Miniconda2-4.6.14-Linux-x86_64.sh && bash Miniconda2-4.6.14-Linux-x86_64.sh -p /miniconda2 -b
ENV PATH=/miniconda2/bin:$PATH
RUN conda update -y conda && rm Miniconda2-4.6.14-Linux-x86_64.sh
RUN conda install -y -c bioconda cutadapt=1.18
RUN conda install -y -c bioconda seqtk=1.3
RUN conda install -y -c bioconda picard=2.20.2

ENV PATH=/usr/local/bin/bwa/:$PATH
ENV PATH=/usr/local/bin/samtools/:$PATH
ENV PATH=/usr/local/bin/samstat/src/:$PATH

ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
RUN mkdir -p scripts
COPY scripts/*.sh scripts/
COPY scripts/BrdUIP-analysis/*.sh scripts/
RUN mkdir -p util
COPY util/ util/

ENV PATH=/usr/local/bin/scripts/:$PATH

RUN chmod +x scripts/*.sh
RUN chmod +x util/*.sh

CMD ["pipeline-list.sh"]

