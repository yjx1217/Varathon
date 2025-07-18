#!/bin/bash
# last update: 2024.05.24
set -e -o pipefail

#########################
VARATHON_HOME=$(pwd)
BUILD="build"
mainland_china_installation="no";
#########################

timestamp () {
  date +"%F %T"
}

clean () {
    dir=$1
    if [ -d $dir ] 
    then
	echo "remove previously failed installation in $BUILD/$dir"
	rm -rf $dir
    fi
}
	
#download () {
#  url=$1
#  download_location=$2
#  echo "Downloading $url to $download_location"
#  wget -c --no-check-certificate $url -O $download_location
#}
download () {
  url=$1
  download_location=$2
  echo "Downloading $url to $download_location"
  if [[ -f $1 ]];then
    # if the tool has been downloaded and deposited in the local, then copy that to destination and install it
    echo -n "The tool exists in the local: $1"
    cp $1 $2
  else
    wget -c --no-check-certificate --max-redirect=30 $url -O $download_location
  fi
}

clone () {
  url=$1
  dir=$(basename $url)
  echo "run clone for \"git clone $url\""
  git clone $url --depth 1
  cd $dir
  git fetch --unshallow
}

tidy_version () { 
    echo "$1" | awk -F. '{ printf("%d%03d%03d%03d\n", $1,$2,$3,$4); }';
}

check_installed () {
    if [ -e "$1/installed" ]; then
        echo "installed"
    else
        echo ""
    fi
}

note_installed () {
    touch "$1/installed"
}

echo ""
echo ""
echo "##################################################################"
echo "###                                                            ###"
echo "###                  Welcome to Varathon                       ###"
echo "###                                                            ###"
echo "##################################################################"
echo ""
echo ""
echo "[$(timestamp)] Installation starts ..."
echo ""


if [ -z "$MAKE_JOBS" ]
then
    echo "Defaulting to 2 concurrent jobs when executing make. Override with MAKE_JOBS=<NUM>"
    MAKE_JOBS=2
fi

if [ ! -z "$INSTALL_DEPS" ]; then
    echo "Installing Varathon build dependencies for Debian/Ubuntu."
    echo "sudo privileges are required and you will be prompted to enter your password"
    sudo apt-get update
    xargs -a debiandeps sudo apt-get install -y
fi

while getopts ":hc" opt
do
    case "${opt}" in
        h)
	    echo "Usage:"
	    echo "bash install_dependencies.sh"
            echo "When running installation within mainland China, please run this script with the '-c' option >"
	    echo "bash install_dependencies.sh -c"
	    echo "";;
        c)
	    echo "Detected the '-c' option >"
	    echo "Set the conda source to the mirror sites in mainland_china for more stable internet access" 
            mainland_china_installation="yes";;
    esac
done
echo "";

# genome simulation
SIMUG_VERSION="1.0.1" # 
SIMUG_GITHUB_COMMIT_VERSION="16e4eea" # committed on 2022.05.25
SIMUG_DOWNLOAD_URL="https://github.com/yjx1217/simuG.git"

# reads retrieving/simulation/processing
SRA_VERSION="3.0.2" # released on 2022.12.12
SRA_DOWNLOAD_URL="https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/${SRA_VERSION}/sratoolkit.${SRA_VERSION}-centos_linux64.tar.gz"

NANOPLOT_VERSION="1.41.0" # released on 2022.12.31

ART_VERSION="mountrainier2016.06.05" # released on 2016.06.05
ART_DOWNLOAD_URL="https://www.niehs.nih.gov/research/resources/assets/docs/artbin${ART_VERSION}linux64.tgz"

PBSIM3_VERSION="f3e5f1c" # committed on 2022.08.01
PBSIM3_DOWNLOAD_URL="https://github.com/yukiteruono/pbsim3.git"

TRIMMOMATIC_VERSION="0.39" # released on 
TRIMMOMATIC_DOWNLOAD_URL="http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-${TRIMMOMATIC_VERSION}.zip"


PORECHOP_VERSION="0.2.4" # 
PORECHOP_GITHUB_COMMIT_VERSION="109e437" # committed on 2018.10.19
PORECHOP_DOWNLOAD_URL="https://github.com/rrwick/Porechop.git"

FILTLONG_VERSION="0.2.1" #
FILTLONG_GITHUB_COMMIT_VERSION="c56f02b" # committed on 2021.07.30
FILTLONG_DOWNLOAD_URL="https://github.com/rrwick/Filtlong.git"

# alignment building/processing
BWA_VERSION="0.7.17" # released on 2017.10.23
BWA_GITHUB_COMMIT_VERSION="139f68f" # committed on 2021.07.30                                                                              
BWA_DOWNLOAD_URL="https://github.com/lh3/bwa.git"
#BWA_DOWNLOAD_URL="https://github.com/lh3/bwa/releases/download/v${BWA_VERSION}/bwa-${BWA_VERSION}.tar.bz2"

BWAMEM2_VERSION="2.2.1" # released on 2021.03.17
BWAMEM2_DOWNLOAD_URL="https://github.com/bwa-mem2/bwa-mem2/releases/download/v${BWAMEM2_VERSION}/bwa-mem2-${BWAMEM2_VERSION}_x64-linux.tar.bz2"

LAST_VERSION="1452" # released on 2023.01.30
LAST_DOWNLOAD_URL="https://gitlab.com/mcfrith/last/-/archive/${LAST_VERSION}/last-${LAST_VERSION}.zip"

NGMLR_VERSION="0.2.7" # released on 2018.06.25
NGMLR_DOWNLOAD_URL="https://github.com/philres/ngmlr/releases/download/v${NGMLR_VERSION}/ngmlr-${NGMLR_VERSION}-linux-x86_64.tar.gz"

MINIMAP2_VERSION="2.24" # released on 2021.12.27
MINIMAP2_GITHUB_COMMIT_VERSION="1d3c3ee" # committed on 2023.02.14
MINIMAP2_DOWNLOAD_URL="https://github.com/lh3/minimap2/releases/download/v${MINIMAP2_VERSION}/minimap2-${MINIMAP2_VERSION}_x64-linux.tar.bz2"

WINNOWMAP_VERSION="2.03" # released on 2021.12.26
WINNOWMAP_GITHUB_COMMIT_VERSION="45078cf" # committed on 2022.04.05
WINNOWMAP_DOWNLOAD_URL="https://github.com/marbl/Winnowmap/archive/refs/tags/v${WINNOWMAP_VERSION}.tar.gz"

LRA_VERSION="1.3.4" # released on 2022.12.14
LRA_DOWNLOAD_URL="https://github.com/ChaissonLab/LRA/archive/refs/tags/${LRA_VERSION}.tar.gz"

PBMM2_VERSION="1.10.0" # released on 2023.01.10
# distributed via conda

GRAPHMAP_VERSION="0.5.2"
GRAPHMAP_GITHUB_COMMIT_VERSION="eb8c75d"
GRAPHMAP_DOWNLOAD_URL="https://github.com/isovic/graphmap"

#GRAPHMAP2_VERSION="0.6.4"
#GRAPHMAP2_DOWNLOAD_URL="https://github.com/lbcb-sci/graphmap2/releases/download/v0.6.4/prebuild.zip"

SAMTOOLS_VERSION="1.17" # released on 2023.02.21
SAMTOOLS_DOWNLOAD_URL="https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2"

HTSLIB_VERSION="1.17" # released on 2023.02.21
HTSLIB_DOWNLOAD_URL="https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2"

#SAMBAMBA_VERSION="1.14" # released on 2021.10.22
#SAMBAMBA_DOWNLOAD_URL="https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2"

PICARD_VERSION="2.27.5" # released on 2022.10.07
PICARD_DOWNLOAD_URL="https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard.jar"

#GATK3_VERSION="3.6-6-g965413b" # 
#GATK3_DOWLOAD_URL="https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/${SRA_VERSION}/GenomeAnalysisTK.jar"
GATK3_DOWNLOAD_URL="https://github.com/yjx1217/GATK3_Archive.git"

GEMTOOLS_VERSION="1.7.1"
GEMTOOLS_DOWNLOAD_URL="http://barnaserver.com/gemtools/releases/GEMTools-static-i3-${GEMTOOLS_VERSION}.tar.gz"

# variant calling
GATK4_VERSION="4.3.0.0" # released on 2022.10.13
GATK4_DOWNLOAD_URL="https://github.com/broadinstitute/gatk/releases/download/${GATK4_VERSION}/gatk-${GATK4_VERSION}.zip"

DEEPVARIANT_VERSION="1.4.0" # released on 2023.03.01

FREEBAYES_VERSION="1.3.6" # released on 2023.03.11
#FREEBAYES_SOURCE_DOWNLOAD_URL="https://github.com/freebayes/freebayes/archive/refs/tags/v${FREEBAYES_VERSION}.tar.gz"

XATLAS_VERSION="0.3"
XATLAS_GITHUB_COMMIT_VERSION="7324754"
XATLAS_DOWNLOAD_URL="https://github.com/jfarek/xatlas/archive/refs/tags/v${XATLAS_VERSION}.tar.gz"

VCFLIB_VERSION="1.0.3" # released on 2023.04.26

CLAIR3_VERSION="1.0.1" # released on 2023.04.26
CLAIR3ILLUMINA_VERSION="1.0.1" # released on 2023.03.08

LONGSHOT_VERSION="0.4.5" # released on 2022.05.03

#WHATSHAP_VERSION="1.2.1" # released on 2021.12.08

CNVKIT_VERSION="0.9.9" # released on 2021.05.29

FREEC_VERSION="11.6" # released on 2020.05.29
#FREEC_VERSION="11.4" # released on 2020.05.29
FREEC_DOWNLOAD_URL="https://github.com/BoevaLab/FREEC/archive/v${FREEC_VERSION}.tar.gz"

MANTA_VERSION="1.6.0" # released on 2019.07.10
MANTA_DOWNLOAD_URL="https://github.com/Illumina/manta/releases/download/v${MANTA_VERSION}/manta-${MANTA_VERSION}.centos6_x86_64.tar.bz2"

SVABA_VERSION="1.1.0" # released on 2022.03.03


DELLY_VERSION="1.1.6" # released on 2022.11.08
DELLY_DOWNLOAD_URL="https://github.com/dellytools/delly/releases/download/v${DELLY_VERSION}/delly_v${DELLY_VERSION}_linux_x86_64bit"

SVIM_VERSION="2.0.0" # released on 2021.06.18
SVIM_GITHUB_COMMIT_VERSION="a9b0985" # committed on 2021.06.29

SNIFFLES_VERSION="2.0.7" # released on 2022.07.25
SNIFFLES_DOWNLOAD_URL="https://github.com/fritzsedlazeck/Sniffles/archive/refs/tags/v${SNIFFLES_VERSION}.tar.gz"

PBSV_VERSION="2.8.0" # released on 2021.06.02
PBSV_GITHUB_COMMIT_VERSION="b9c5504"
PBSV_DOWNLOAD_URL="https://github.com/PacificBiosciences/pbsv"

PICKY_VERSION="0.2.a" # released on 2018.07.17   
PICKY_GITHUB_COMMIT_VERSION="34b85ac" # committed on 2018.07.17
PICKY_DOWNLOAD_URL="https://github.com/TheJacksonLaboratory/Picky"

NANOSV_VERSION="1.2.4" # released on 2019.04.09
NANOSV_GITHUB_COMMIT_VERSION="c1ae30c" # committed on 2019.04.09
NANOSV_DOWNLOAD_URL="https://github.com/mroosmalen/nanosv"

NANOVAR_VERSION="1.4.1" # released on 2022.01.17

CUTESV_VERSION="2.0.2" # released on 2022.11.02

DEBREAK_VERSION="1.3" # released on 2022.09.02
DEBREAK_GITHUB_COMMIT_VERSION="56b3f94"
DEBREAK_DOWNLOAD_URL="https://github.com/ChongLab/DeBreak"

DUPHOLD_VERSION="0.2.3"
DUPHOLD_GITHUB_COMMIT_VERSION="e14d6ba"
DUPHOLD_DOWNLOAD_URL="https://github.com/brentp/duphold/releases/download/v0.2.3/duphold"

JASMINESV_VERSION="1.1.5" # released on 2022.04.28

USABLEVCF_GITHUB_COMMIT_VERSION="07c39ae" # commited on 2020.05.16

TRUVARI_VERSION="4.0.0" # released on 2023.03.14

SANSA_VERSION="0.0.8"
SANSA_DOWNLOAD_URL="https://github.com/dellytools/sansa/releases/download/v${SANSA_VERSION}/sansa_v${SANSA_VERSION}_linux_x86_64bit"

SAMPLOT_VERSION="1.3.0" # released on 2021.10.14


# variant processing
VT_VERSION="" # we use the github commit version below
VT_GITHUB_COMMIT_VERSION="f6d2b5d" # "ab6c13c" # committed on 2021.05.04
VT_DOWNLOAD_URL="https://github.com/atks/vt"


BCFTOOLS_VERSION="1.17" # released on 2018.07.18
BCFTOOLS_DOWNLOAD_URL="https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2"

MINICONDA3_VERSION="py39_25.5.1-0" # released on 2025.06.24
if [[ "$mainland_china_installation" == "no" ]]
then
    MINICONDA3_DOWNLOAD_URL="https://repo.anaconda.com/miniconda/Miniconda3-${MINICONDA3_VERSION}-Linux-x86_64.sh"
else
    MINICONDA3_DOWNLOAD_URL="https://mirrors.tuna.tsinghua.edu.cn/anaconda/miniconda/Miniconda3-${MINICONDA3_VERSION}-Linux-x86_64.sh"
fi

PYTHON3_VERSION="3.8"

BEDTOOLS_VERSION="2.31.1" # released on 2023.11.08
BEDTOOLS_DOWNLOAD_URL="https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLS_VERSION}/bedtools-${BEDTOOLS_VERSION}.tar.gz"

BLAST_VERSION="2.2.31"
BLAST_DOWNLOAD_URL="https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${BLAST_VERSION}/ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz"

RMBLAST_VERSION="2.2.28"
RMBLAST_DOWNLOAD_URL="https://ftp.ncbi.nlm.nih.gov/blast/executables/rmblast/${RMBLAST_VERSION}/ncbi-rmblastn-${RMBLAST_VERSION}-x64-linux.tar.gz"

VEP_VERSION="109.3"

PARALLEL_VERSION="20180722" # released on 2018.07.22
PARALLEL_DOWNLOAD_URL="http://ftp.gnu.org/gnu/parallel/parallel-${PARALLEL_VERSION}.tar.bz2"

# MUMMER4_VERSION="4.0.0beta2" # released on 2017.10.14
# MUMMER4_DOWNLOAD_URL="https://github.com/gmarcais/mummer/releases/download/v${MUMMER_VERSION}/mummer-${MUMMER_VERSION}.tar.gz"

# GNUPLOT_VERSION="4.6.6"
# GNUPLOT_DOWNLOAD_URL="https://sourceforge.net/projects/gnuplot/files/gnuplot/${GNUPLOT_VERSION}/gnuplot-${GNUPLOT_VERSION}.tar.gz"

if [ -d $BUILD ]
then
    echo ""
    echo "[$(timestamp)] Detected previously generated $BUILD directory."
else
    echo "[$(timestamp)] Create the new $BUILD directory."
    mkdir $BUILD
    echo ""
fi

cd $BUILD
build_dir=$(pwd)

# Downloading all the dependencies
echo ""
echo "[$(timestamp)] Download and install all the dependencies ..."
echo ""

# ---------- set Perl & R environment variables -------------

PYTHONPATH=""
PERL5LIB="$build_dir:$PERL5LIB"
PERL5LIB="$build_dir/cpanm/perlmods/lib/perl5:$PERL5LIB"
R_LIBS="$build_dir/R_libs:$R_LIBS"

echo "[$(timestamp)] Installing Perl modules ..."
cpanm_dir="$build_dir/cpanm"
if [ -z $(check_installed $cpanm_dir) ]; then
    mkdir -p  $cpanm_dir
    cd $cpanm_dir
    # wget -c --no-check-certificate -O - https://cpanmin.us/ > cpanm
    # work around for the unstable downloading issue
    cp $VARATHON_HOME/misc/cpanm .
    chmod +x cpanm
    mkdir -p perlmods
    $cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed Test::More@1.302086
    $cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed Env@1.04
    $cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed Statistics::Descriptive@3.0612
    $cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed Statistics::Descriptive::Discrete@0.07
    $cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed Math::Random@0.72
    $cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed Math::Round@0.07
    note_installed $cpanm_dir
fi    

echo "[$(timestamp)] Installing R libraries ..."
rlib_dir="$build_dir/R_libs"
mkdir -p $rlib_dir
cd $rlib_dir
R_VERSION=$(R --version |head -1 |cut -d " " -f 3)

if [ -z $(check_installed "$rlib_dir/optparse") ]; then
    clean "$rlib_dir/optparse"
    R -e "install.packages(\"optparse\", repos=\"http://cran.rstudio.com/\", lib=\"$build_dir/R_libs/\")"
    note_installed "$rlib_dir/optparse"
fi

if [ -z $(check_installed "$rlib_dir/ggplot2") ]; then
    clean "$rlib_dir/ggplot2"
    R -e "install.packages(\"ggplot2\", repos=\"http://cran.rstudio.com/\", lib=\"$build_dir/R_libs/\")"
    note_installed "$rlib_dir/ggplot2"
fi

if [ -z $(check_installed "$rlib_dir/scales") ]; then
    clean "$rlib_dir/scales"
    R -e "install.packages(\"scales\", repos=\"http://cran.rstudio.com/\", lib=\"$build_dir/R_libs/\")"
    note_installed "$rlib_dir/scales"
fi

if [ -z $(check_installed "$rlib_dir/viridis") ]; then
    clean "$rlib_dir/viridis"
    R -e "install.packages(\"viridis\", repos=\"http://cran.rstudio.com/\", lib=\"$build_dir/R_libs/\")"
    note_installed "$rlib_dir/viridis"
fi

if [ -z $(check_installed "$rlib_dir/DNAcopy") ]; then
    clean "$rlib_dir/DNAcopy"
    if [ $(tidy_version "$R_VERSION") -ge $(tidy_version "3.6.0") ]
    then
	echo "R_VERSION=$R_VERSION, use the new bioconductor installation protocol"
	R -e ".libPaths(\"$build_dir/R_libs/\");install.packages(\"BiocManager\", repos=\"http://cran.rstudio.com/\", lib=\"$build_dir/R_libs/\");BiocManager::install(\"DNAcopy\", lib=\"$build_dir/R_libs/\")"
    else
	echo "R_VERSION=$R_VERSION, use the old bioconductor installation protocol"
	R -e ".libPaths(\"$build_dir/R_libs/\");source(\"https://bioconductor.org/biocLite.R\");biocLite(\"DNAcopy\", lib=\"$build_dir/R_libs/\", type = \"source\")"
    fi
    note_installed "$rlib_dir/DNAcopy"
fi


# install dependencies

# # ------------- Miniconda2 --------------------
# miniconda2_dir="$build_dir/miniconda2/bin"
# if [ -z $(check_installed $miniconda2_dir) ]; then
#     cd $build_dir
#     download $MINICONDA2_DOWNLOAD_URL "Miniconda2-${MINICONDA2_VERSION}-Linux-x86_64.sh"
#     bash Miniconda2-${MINICONDA2_VERSION}-Linux-x86_64.sh -b -p $build_dir/miniconda2
#     if [[ "$mainland_china_installation" == "yes" ]]
#     then
# 	$miniconda2_dir/conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/pkgs/main
# 	$miniconda2_dir/conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/pkgs/free
# 	$miniconda2_dir/conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/pkgs/pro
# 	$miniconda2_dir/conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/pkgs/msys2
# 	$miniconda2_dir/conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/cloud/conda-forge
# 	$miniconda2_dir/conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/cloud/bioconda

# 	$miniconda2_dir/conda config --add channels https://anaconda.mirrors.sjtug.sjtu.edu.cn/pkgs/main
# 	$miniconda2_dir/conda config --add channels https://anaconda.mirrors.sjtug.sjtu.edu.cn/pkgs/free
# 	$miniconda2_dir/conda config --add channels https://anaconda.mirrors.sjtug.sjtu.edu.cn/pkgs/pro
# 	$miniconda2_dir/conda config --add channels https://anaconda.mirrors.sjtug.sjtu.edu.cn/pkgs/msys2
# 	$miniconda2_dir/conda config --add channels https://anaconda.mirrors.sjtug.sjtu.edu.cn/cloud/conda-forge
# 	$miniconda2_dir/conda config --add channels https://anaconda.mirrors.sjtug.sjtu.edu.cn/cloud/bioconda
#     else 
# 	$miniconda2_dir/conda config --add channels defaults
# 	$miniconda2_dir/conda config --add channels bioconda
# 	$miniconda2_dir/conda config --add channels conda-forge
#     fi
#     $miniconda2_dir/conda config --set show_channel_urls yes
#     $miniconda2_dir/conda config --set channel_priority flexible
#     rm Miniconda2-${MINICONDA2_VERSION}-Linux-x86_64.sh 
#     note_installed $miniconda2_dir
# else
#     if [[ "$mainland_china_installation" == "yes" ]]
#     then
# 	$miniconda2_dir/conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/pkgs/main
# 	$miniconda2_dir/conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/pkgs/free
# 	$miniconda2_dir/conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/pkgs/pro
# 	$miniconda2_dir/conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/pkgs/msys2
# 	$miniconda2_dir/conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/cloud/conda-forge
# 	$miniconda2_dir/conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/cloud/bioconda

# 	$miniconda2_dir/conda config --add channels https://anaconda.mirrors.sjtug.sjtu.edu.cn/pkgs/main
# 	$miniconda2_dir/conda config --add channels https://anaconda.mirrors.sjtug.sjtu.edu.cn/pkgs/free
# 	$miniconda2_dir/conda config --add channels https://anaconda.mirrors.sjtug.sjtu.edu.cn/pkgs/pro
# 	$miniconda2_dir/conda config --add channels https://anaconda.mirrors.sjtug.sjtu.edu.cn/pkgs/msys2
# 	$miniconda2_dir/conda config --add channels https://anaconda.mirrors.sjtug.sjtu.edu.cn/cloud/conda-forge
# 	$miniconda2_dir/conda config --add channels https://anaconda.mirrors.sjtug.sjtu.edu.cn/cloud/bioconda

#     else 
# 	$miniconda2_dir/conda config --add channels defaults
# 	$miniconda2_dir/conda config --add channels bioconda
# 	$miniconda2_dir/conda config --add channels conda-forge
#     fi
#     $miniconda2_dir/conda config --set show_channel_urls yes
#     $miniconda2_dir/conda config --set channel_priority flexible
# fi


# ------------- Miniconda3 --------------------
echo "[$(timestamp)] Installing miniconda3 ..."
miniconda3_dir="$build_dir/miniconda3/bin"
if [ -z $(check_installed $miniconda3_dir) ]; then
    cd $build_dir
    clean "$build_dir/miniconda3"
    download $MINICONDA3_DOWNLOAD_URL "Miniconda3-${MINICONDA3_VERSION}-Linux-x86_64.sh"
    bash Miniconda3-${MINICONDA3_VERSION}-Linux-x86_64.sh -b -p $build_dir/miniconda3
    if [[ "$mainland_china_installation" == "yes" ]]
    then
	$miniconda3_dir/conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/pkgs/main
	$miniconda3_dir/conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/pkgs/free
	$miniconda3_dir/conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/pkgs/pro
	$miniconda3_dir/conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/pkgs/msys2
	$miniconda3_dir/conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/cloud/bioconda
	$miniconda3_dir/conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/cloud/conda-forge

	$miniconda3_dir/conda config --add channels https://anaconda.mirrors.sjtug.sjtu.edu.cn/pkgs/main
	$miniconda3_dir/conda config --add channels https://anaconda.mirrors.sjtug.sjtu.edu.cn/pkgs/free
	$miniconda3_dir/conda config --add channels https://anaconda.mirrors.sjtug.sjtu.edu.cn/pkgs/pro
	$miniconda3_dir/conda config --add channels https://anaconda.mirrors.sjtug.sjtu.edu.cn/pkgs/msys2
	$miniconda3_dir/conda config --add channels https://anaconda.mirrors.sjtug.sjtu.edu.cn/cloud/bioconda
	$miniconda3_dir/conda config --add channels https://anaconda.mirrors.sjtug.sjtu.edu.cn/cloud/conda-forge
    else 
        $miniconda3_dir/conda config --add channels defaults
        $miniconda3_dir/conda config --add channels bioconda
        $miniconda3_dir/conda config --add channels conda-forge
    fi
    $miniconda3_dir/conda config --set ssl_verify false
    $miniconda3_dir/conda config --set show_channel_urls yes
    $miniconda3_dir/conda config --set channel_priority flexible
    $miniconda3_dir/conda config --set report_errors false
    $miniconda3_dir/conda clean -i 
    # $miniconda3_dir/conda install pip numpy scipy cython pysam 
    cd $build_dir
    rm Miniconda3-${MINICONDA3_VERSION}-Linux-x86_64.sh 
    note_installed $miniconda3_dir
else
    if [[ "$mainland_china_installation" == "yes" ]]
    then
	$miniconda3_dir/conda clean -i 

	$miniconda3_dir/conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/pkgs/main
	$miniconda3_dir/conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/pkgs/free
	$miniconda3_dir/conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/pkgs/pro
	$miniconda3_dir/conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/pkgs/msys2
	$miniconda3_dir/conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/cloud/bioconda
	$miniconda3_dir/conda config --add channels https://mirrors.bfsu.edu.cn/anaconda/cloud/conda-forge

	$miniconda3_dir/conda config --add channels https://anaconda.mirrors.sjtug.sjtu.edu.cn/pkgs/main
	$miniconda3_dir/conda config --add channels https://anaconda.mirrors.sjtug.sjtu.edu.cn/pkgs/free
	$miniconda3_dir/conda config --add channels https://anaconda.mirrors.sjtug.sjtu.edu.cn/pkgs/pro
	$miniconda3_dir/conda config --add channels https://anaconda.mirrors.sjtug.sjtu.edu.cn/pkgs/msys2
	$miniconda3_dir/conda config --add channels https://anaconda.mirrors.sjtug.sjtu.edu.cn/cloud/bioconda
	$miniconda3_dir/conda config --add channels https://anaconda.mirrors.sjtug.sjtu.edu.cn/cloud/conda-forge
    else 
	$miniconda3_dir/conda config --add channels defaults
	$miniconda3_dir/conda config --add channels bioconda
	$miniconda3_dir/conda config --add channels conda-forge
    fi
    $miniconda3_dir/conda config --set ssl_verify false
    $miniconda3_dir/conda config --set show_channel_urls yes
    $miniconda3_dir/conda config --set channel_priority flexible
    $miniconda3_dir/conda config --set report_errors false
    $miniconda3_dir/conda clean -i 
fi


# ------------- simuG -------------------
echo ""
echo "[$(timestamp)] Installing simuG ..."
simuG_dir="$build_dir/simuG"
if [ -z $(check_installed $simuG_dir) ]; then
    cd $build_dir
    clean "$build_dir/simuG"
    echo "Download simuG-v${SIMUG_VERSION}"
    git clone $SIMUG_DOWNLOAD_URL
    cd simuG
    git checkout -f -q $SIMUG_GITHUB_COMMIT_VERSION
    note_installed $simuG_dir
fi

# ------------- SRA Toolkit -------------------
echo ""
echo "[$(timestamp)] Installing SRA toolkit ..."
sra_dir="$build_dir/sratoolkit.${SRA_VERSION}-centos_linux64/bin"
if [ -z $(check_installed $sra_dir) ]; then
    cd $build_dir
    clean "$build_dir/sratoolkit.${SRA_VERSION}-centos_linux64"
    echo "Download SRAtoolkit-v${SRA_VERSION}"
    download $SRA_DOWNLOAD_URL sratoolkit.${SRA_VERSION}-centos_linux64.tar.gz
    tar -xzf sratoolkit.${SRA_VERSION}-centos_linux64.tar.gz
    rm sratoolkit.${SRA_VERSION}-centos_linux64.tar.gz
    note_installed $sra_dir
fi

# --------------- Nanoplot --------------------
echo ""
echo "[$(timestamp)] Installing nanoplot ..."
nanoplot_dir="$build_dir/nanoplot_conda_env/bin"
if [ -z $(check_installed $nanoplot_dir) ]; then
    cd $build_dir
    clean "$build_dir/nanoplot_conda_env"
    $miniconda3_dir/conda create -y -p $build_dir/nanoplot_conda_env python=${PYTHON3_VERSION}
    source $miniconda3_dir/activate $build_dir/nanoplot_conda_env
    # if [[ "$mainland_china_installation" == "yes" ]]; then
    # 	$miniconda3_dir/conda install -y nanoplot=${NANOPLOT_VERSION}

    # else
    # 	$miniconda3_dir/conda install -y -c bioconda nanoplot=${NANOPLOT_VERSION}
    # fi
    $miniconda3_dir/pip3 install NanoPlot==${NANOPLOT_VERSION}
    source $miniconda3_dir/deactivate
    note_installed $nanoplot_dir
fi

# ------------- ART -------------------
echo ""
echo "[$(timestamp)] Installing ART ..."
art_dir="$build_dir/art_bin_MountRainier"
if [ -z $(check_installed $art_dir) ]; then
    cd $build_dir
    clean "$build_dir/art_bin_MountRainier"
    echo "Download ART-v${ART_VERSION}"
    download $ART_DOWNLOAD_URL artbin${ART_VERSION}linux64.tgz
    tar -zxf artbin${ART_VERSION}linux64.tgz
    rm artbin${ART_VERSION}linux64.tgz
    note_installed $art_dir
fi

# ------------- PBSIM3 -------------------
echo ""
echo "[$(timestamp)] Installing pbsim3 ..."
pbsim3_dir="$build_dir/pbsim3/build"
if [ -z $(check_installed $pbsim3_dir) ]; then
    cd $build_dir
    clean "$build_dir/pbsim3"
    echo "Download pbsim3-v${PBSIM3_VERSION}"
    git clone $PBSIM3_DOWNLOAD_URL
    cd $build_dir/pbsim3
    git checkout -f -q $PBSIM3_GITHUB_COMMIT_VERSION
    mkdir build
    ./configure --prefix="$(pwd)/build"
    make
    make install
    note_installed $pbsim3_dir
fi

# # # -------------- SimLoRD ----------------
# echo ""
# echo "[$(timestamp)] Installing SimLoRd ..."
# simlord_dir="$build_dir/simlord_conda_env/bin"
# if [ -z $(check_installed $simlord_dir) ]; then
#     cd $build_dir
#     clean "$build_dir/simlord_conda_env"
#     echo "Download simlord-v${SIMLORD_VERSION}"
#     $miniconda3_dir/conda create -y -p $build_dir/simlord_conda_env python=${PYTHON3_VERSION} pip
#     source $miniconda3_dir/activate $build_dir/simlord_conda_env
#     # if [[ "$mainland_china_installation" == "yes" ]]; then
#     # 	$miniconda3_dir/conda install -y simlord=${SIMLORD_VERSION}
#     # else
#     # 	$miniconda3_dir/conda install -y -c bioconda simlord=${SIMLORD_VERSION}
#     # fi
#     $miniconda3_dir/pip3 install numpy==1.21.4
#     $miniconda3_dir/pip3 install cython==0.29.24
#     $miniconda3_dir/pip3 install scipy==1.7.2
#     $miniconda3_dir/pip3 install pysam==0.17.0
#     $miniconda3_dir/pip3 install dinopy==2.2.0
#     $miniconda3_dir/pip3 install "simlord==${SIMLORD_VERSION}"
#     source $miniconda3_dir/deactivate
#     note_installed $simlord_dir
# fi

# # -------------- NanoSim ----------------
# echo ""
# echo "[$(timestamp)] Installing NanoSim ..."
# nanosim_dir="$build_dir/nanosim_conda_env/bin"
# if [ -z $(check_installed $nanosim_dir) ]; then
#     cd $build_dir
#     clean "$build_dir/nanosim_conda_env"
#     echo "Download NanoSim-v${NANOSIM_VERSION}"
#     $miniconda3_dir/conda create -y -p $build_dir/nanosim_conda_env python=${PYTHON3_VERSION} # pip numpy scipy cython
#     source $miniconda3_dir/activate $build_dir/nanosim_conda_env
#     if [[ "$mainland_china_installation" == "yes" ]]; then
# 	$miniconda3_dir/conda install -y nanosim=${NANOSIM_VERSION}
#     else
# 	$miniconda3_dir/conda install -y -c bioconda nanosim=${NANOSIM_VERSION}

#     fi
#     source $miniconda3_dir/deactivate
#     note_installed $nanosim_dir
# fi

# # -------------- DeepSimulator ----------------
# deepsimulator_dir="$build_dir/DeepSimulator"
# if [ -z $(check_installed $deepsimulator_dir) ]; then
#     cd $build_dir
#     clean $deepsimulator_dir
#     echo "Download DeepSimulator-v${DEEPSIMULATOR_GITHUB_COMMIT_VERSION}"
#     # git clone $DEEPSIMULATOR_DOWNLOAD_URL
#     clone $DEEPSIMULATOR_DOWNLOAD_URL
#     cd $deepsimulator_dir
#     git checkout -f -q $DEEPSIMULATOR_GITHUB_COMMIT_VERSION
#     $miniconda2_dir/conda remove --name tensorflow_cdpm --all -y
#     $miniconda2_dir/conda create --name tensorflow_cdpm python=2.7 -y
#     source $miniconda2_dir/activate tensorflow_cdpm
#     $miniconda2_dir/conda install -y -c anaconda scikit-learn=0.20.3
#     $miniconda2_dir/pip install numpy==1.13.1
#     if [[ "$mainland_china_installation" == "yes" ]]
#     then
# 	$miniconda2_dir/pip install tensorflow==1.2.1 -i https://pypi.tuna.tsinghua.edu.cn/simple/
#     else
# 	$miniconda2_dir/pip install tensorflow==1.2.1
#     fi
#     $miniconda2_dir/pip install tflearn==0.3.2
#     $miniconda2_dir/pip install tqdm==4.19.4
#     $miniconda2_dir/pip install scipy==0.18.1
#     $miniconda2_dir/pip install h5py==2.7.1
#     $miniconda2_dir/pip install biopython==1.74
#     source $miniconda2_dir/deactivate
#     $miniconda2_dir/conda remove --name basecall --all -y
#     $miniconda2_dir/conda create --name basecall python=${PYTHON3_VERSION} -y
#     source $miniconda2_dir/activate basecall
#     #-> 2. install basecaller
#     #--| 2.1 install albacore_2.3.1
#     cd base_caller/albacore_2.3.1/
#     wget -q https://mirror.oxfordnanoportal.com/software/analysis/ont_albacore-2.3.1-cp36-cp36m-manylinux1_x86_64.whl
#     pip install ont_albacore-2.3.1-cp36-cp36m-manylinux1_x86_64.whl
#     rm -f ont_albacore-2.3.1-cp36-cp36m-manylinux1_x86_64.whl
#     cd ../../
#     #--| 2.2 install guppy_3.1.5
#     cd base_caller/guppy_3.1.5/
#     wget -q https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy-cpu_3.1.5_linux64.tar.gz
#     tar xzf ont-guppy-cpu_3.1.5_linux64.tar.gz
#     rm -f ont-guppy-cpu_3.1.5_linux64.tar.gz
#     cd ../../
#     source $miniconda2_dir/deactivate
#     cd $deepsimulator_dir
#     note_installed $deepsimulator_dir
# fi

# --------------- Trimmomatic -----------------
echo ""
echo "[$(timestamp)] Installing trimmomatic ..."
trimmomatic_dir="$build_dir/Trimmomatic-${TRIMMOMATIC_VERSION}"
if [ -z $(check_installed $trimmomatic_dir) ]; then
    cd $build_dir
    clean "$build_dir/Trimmomatic-${TRIMMOMATIC_VERSION}"
    echo "Download Trimmomatic-v${TRIMMOMATIC_VERSION}"
    download $TRIMMOMATIC_DOWNLOAD_URL "Trimmomatic-${TRIMMOMATIC_VERSION}.zip"
    unzip Trimmomatic-${TRIMMOMATIC_VERSION}.zip
    cd $trimmomatic_dir
    chmod 755 trimmomatic-${TRIMMOMATIC_VERSION}.jar
    ln -s trimmomatic-${TRIMMOMATIC_VERSION}.jar trimmomatic.jar 
    cat $trimmomatic_dir/adapters/*.fa |sed 's/>/\n>/' > adapters.fa
    cd $build_dir
    rm Trimmomatic-${TRIMMOMATIC_VERSION}.zip
    note_installed $trimmomatic_dir
fi

# --------------- Porechop ------------------
echo ""
echo "[$(timestamp)] Installing Porechop ..."
porechop_dir="$build_dir/porechop_conda_env/bin"
if [ -z $(check_installed $porechop_dir) ]; then
    cd $build_dir
    clean "$build_dir/porechop_conda_env"
    echo "Download Porechop-v${PORECHOP_VERSION}"
    $miniconda3_dir/conda create -y -p $build_dir/porechop_conda_env python=${PYTHON3_VERSION}
    source $miniconda3_dir/activate $build_dir/porechop_conda_env
    if [[ "$mainland_china_installation" == "yes" ]]; then
	$miniconda3_dir/conda install -y porechop=${PORECHOP_VERSION} 
    else
	$miniconda3_dir/conda install -y -c bioconda porechop=${PORECHOP_VERSION} 
    fi
    source $miniconda3_dir/deactivate
    note_installed $porechop_dir
fi

# --------------- Filtlong ------------------
echo ""
echo "[$(timestamp)] Installing Filtlong ..."
filtlong_dir="$build_dir/Filtlong/bin"
if [ -z $(check_installed $filtlong_dir) ]; then
    cd $build_dir
    clean "$build_dir/Filtlong"
    echo "Download Filtlong-v${FILTLONG_VERSION}"
    git clone $FILTLONG_DOWNLOAD_URL
    cd Filtlong
    git checkout -f -q $FILTLONG_GITHUB_COMMIT_VERSION
    make -j $MAKE_JOBS
    note_installed $filtlong_dir
fi

# ------------- BWA -------------------
echo ""
echo "[$(timestamp)] Installing bwa ..."
bwa_dir="$build_dir/bwa"
if [ -z $(check_installed $bwa_dir) ]; then
    cd $build_dir
    clean "$build_dir/bwa"
    echo "Download BWA-v${BWA_VERSION}"
    #download $BWA_DOWNLOAD_URL "bwa-${BWA_VERSION}.tar.bz2"
    #tar -xjf bwa-${BWA_VERSION}.tar.bz2
    git clone $BWA_DOWNLOAD_URL
    cd $bwa_dir
    git checkout -f -q $BWA_GITHUB_COMMIT_VERSION
    make -j $MAKE_JOBS
    cd $build_dir
    #rm bwa-${BWA_VERSION}.tar.bz2
    note_installed $bwa_dir
fi

# ------------- BWA-MEM2 -------------------
echo ""
echo "[$(timestamp)] Installing bwa-mem2 ..."
bwamem2_dir="$build_dir/bwa-mem2-${BWAMEM2_VERSION}_x64-linux"
if [ -z $(check_installed $bwamem2_dir) ]; then
    cd $build_dir
    clean "$build_dir/bwa-mem2-${BWAMEM2_VERSION}_x64-linux"
    echo "Download BWAMEM2-v${BWAMEM2_VERSION}"
    download $BWAMEM2_DOWNLOAD_URL "bwa-mem2-${BWAMEM2_VERSION}.tar.bz2"
    tar -xjf bwa-mem2-${BWAMEM2_VERSION}.tar.bz2
    cd $build_dir
    rm bwa-mem2-${BWAMEM2_VERSION}.tar.bz2
    note_installed $bwamem2_dir
fi

# ------------- LAST -------------------
echo ""
echo "[$(timestamp)] Installing LAST ..." 
last_dir="$build_dir/last-${LAST_VERSION}/bin"
if [ -z $(check_installed $last_dir) ]; then
    cd $build_dir
    if [ -f last-${LAST_VERSION}.zip ]; then
	rm last-${LAST_VERSION}.zip
    fi
    clean "$build_dir/last-${LAST_VERSION}"
    echo "Download LAST-v${LAST_VERSION}"
    download $LAST_DOWNLOAD_URL "last-${LAST_VERSION}.zip"
    unzip "last-${LAST_VERSION}.zip"
    cd "$build_dir/last-${LAST_VERSION}/src"
    make -j $MAKE_JOBS
    cd ..
    cd $build_dir
    rm last-${LAST_VERSION}.zip
    note_installed $last_dir
fi

# ------------- NGMLR -------------------
echo ""
echo "[$(timestamp)] Installing NGMLR ..." 
ngmlr_dir="$build_dir/ngmlr-${NGMLR_VERSION}"
if [ -z $(check_installed $ngmlr_dir) ]; then
    cd $build_dir
    clean "$build_dir/ngmlr-${NGMLR_VERSION}"
    echo "Download NGMLR-v${NGMLR_VERSION}"
    download $NGMLR_DOWNLOAD_URL "ngmlr-${NGMLR_VERSION}.tar.gz"
    tar xvzf "ngmlr-${NGMLR_VERSION}.tar.gz"
    rm ngmlr-${NGMLR_VERSION}.tar.gz
    note_installed $ngmlr_dir
fi

# --------------- minimap2 ------------------
echo ""
echo "[$(timestamp)] Installing minimap2 ..."
minimap2_dir="$build_dir/minimap2-${MINIMAP2_VERSION}_x64-linux"
if [ -z $(check_installed $minimap2_dir) ]; then
    cd $build_dir
    clean "$build_dir/minimap2-${MINIMAP2_VERSION}_x64-linux"
    echo "Download minimap2-v${MINIMAP2_VERSION}"
    download $MINIMAP2_DOWNLOAD_URL "minimap2-${MINIMAP2_VERSION}.tar.bz2"
    tar xvjf minimap2-${MINIMAP2_VERSION}.tar.bz2
    rm minimap2-${MINIMAP2_VERSION}.tar.bz2
    note_installed $minimap2_dir
fi

# --------------- winnowmap ------------------
echo ""
echo "[$(timestamp)] Installing winnowmap ..."
winnowmap_dir="$build_dir/Winnowmap-${WINNOWMAP_VERSION}/bin"
if [ -z $(check_installed $winnowmap_dir) ]; then
    cd $build_dir
    clean "$build_dir/Winnowmap-${WINNOWMAP_VERSION}"
    echo "Download winnowmap-v${WINNOWMAP_VERSION}"
    download $WINNOWMAP_DOWNLOAD_URL "winnowmap-${WINNOWMAP_VERSION}.tar.gz"
    tar xvzf winnowmap-${WINNOWMAP_VERSION}.tar.gz
    cd $build_dir/Winnowmap-${WINNOWMAP_VERSION}
    make -j $MAKE_JOBS
    cd $build_dir
    rm winnowmap-${WINNOWMAP_VERSION}.tar.gz
    note_installed $winnowmap_dir
fi

# --------------- LRA -----------------
echo ""
echo "[$(timestamp)] Installing lra ..."
lra_dir="$build_dir/lra_conda_env/bin"
if [ -z $(check_installed $lra_dir) ]; then
    cd $build_dir
    clean "$build_dir/lra_conda_env"
    echo "Download lra-v${LRA_VERSION}"
    $miniconda3_dir/conda create -y -p $build_dir/lra_conda_env python=${PYTHON3_VERSION}
    source $miniconda3_dir/activate $build_dir/lra_conda_env
    if [[ "$mainland_china_installation" == "yes" ]]; then
	$miniconda3_dir/conda install -y lra=${LRA_VERSION}
    else
	$miniconda3_dir/conda install -y -c bioconda lra=${LRA_VERSION}
    fi
    source $miniconda3_dir/deactivate
    note_installed $lra_dir
fi

# --------------- pbmm2 -----------------
echo ""
echo "[$(timestamp)] Installing pbmm2 ..."
pbmm2_dir="$build_dir/pbmm2_conda_env/bin"
if [ -z $(check_installed $pbmm2_dir) ]; then
    cd $build_dir
    clean "$build_dir/pbmm2_conda_env"
    echo "Download pbmm2-v${PBMM2_VERSION}"
    $miniconda3_dir/conda create -y -p $build_dir/pbmm2_conda_env python=${PYTHON3_VERSION}
    source $miniconda3_dir/activate $build_dir/pbmm2_conda_env
    if [[ "$mainland_china_installation" == "yes" ]]; then
	$miniconda3_dir/conda install -y pbmm2=${PBMM2_VERSION}
    else
	$miniconda3_dir/conda install -y -c bioconda pbmm2=${PBMM2_VERSION}
    fi
    source $miniconda3_dir/deactivate
    note_installed $pbmm2_dir
fi

# --------------- GraphMap ------------------
echo ""
echo "[$(timestamp)] Installing GraphMap ..."
graphmap_dir="$build_dir/graphmap/bin/Linux-x64"
if [ -z $(check_installed $graphmap_dir) ]; then
    cd $build_dir
    clean "$build_dir/graphmap"
    echo "Download graphmap-v${GRAPHMAP_VERSION}"
    git clone $GRAPHMAP_DOWNLOAD_URL
    cd graphmap
    git checkout -f -q $GRAPHMAP_GITHUB_COMMIT_VERSION
    cd codebase
    git clone https://github.com/isovic/seqlib.git
    git clone https://github.com/isovic/argumentparser.git
    git clone https://github.com/isovic/gindex
    # make modules
    cd ..
    make -j $MAKE_JOBS
    note_installed $graphmap_dir
fi

# # --------------- GraphMap2 ------------------
# echo ""
# echo "[$(timestamp)] Installing GraphMap2 ..."
# graphmap2_dir="$build_dir/graphmap2-${GRAPHMAP2_VERSION}"
# if [ -z $(check_installed $graphmap2_dir) ]; then
# cd $build_dir
# clean $graphmap2_dir
# echo "Download graphmap2-v${GRAPHMAP2_VERSION}"
# download $GRAPHMAP2_DOWNLOAD_URL "graphmap2-${GRAPHMAP2_VERSION}.zip"
# unzip graphmap2-${GRAPHMAP2_VERSION}.zip
# mv prebuild graphmap2-${GRAPHMAP2_VERSION}
# cd graphmap2-${GRAPHMAP2_VERSION}
# ln -s graphmap-linux graphmap 
# rm graphmap2-${GRAPHMAP2_VERSION}.zip
# note_installed $graphmap2_dir
# fi

# --------------- samtools -----------------
echo ""
echo "[$(timestamp)] Installing samtools,htslib,tabix ..."
samtools_dir="$build_dir/samtools-${SAMTOOLS_VERSION}"
htslib_dir="$samtools_dir/htslib-${HTSLIB_VERSION}"
tabix_dir="$samtools_dir/htslib-${HTSLIB_VERSION}"
if [ -z $(check_installed $samtools_dir) ]; then
    cd $build_dir
    clean "$build_dir/samtools-${SAMTOOLS_VERSION}"
    echo "Download samtools-v${SAMTOOLS_VERSION}"
    download $SAMTOOLS_DOWNLOAD_URL "samtools-${SAMTOOLS_VERSION}.tar.bz2"
    tar xvjf samtools-${SAMTOOLS_VERSION}.tar.bz2
    cd $samtools_dir
    C_INCLUDE_PATH=""
    ./configure --without-curses;
    make -j $MAKE_JOBS
    cd htslib-${HTSLIB_VERSION}
    ./configure
    make -j $MAKE_JOBS
    cd $build_dir
    rm samtools-${SAMTOOLS_VERSION}.tar.bz2
    note_installed $samtools_dir
fi
PATH="$samtools_dir:$htslib_dir:$tabix_dir:${PATH}"

# --------------- Picard -----------------
echo ""
echo "[$(timestamp)] Installing Picard ..."
picard_dir="$build_dir/Picard-v${PICARD_VERSION}"
if [ -z $(check_installed $picard_dir) ]; then
    cd $build_dir
    clean "$build_dir/Picard-v${PICARD_VERSION}"
    echo "Download Picard-v${PICARD_VERSION}"
    download $PICARD_DOWNLOAD_URL "picard.jar"
    mkdir Picard-v${PICARD_VERSION}
    mv picard.jar $picard_dir
    cd $picard_dir
    chmod 755 picard.jar
    note_installed $picard_dir
fi

# --------------- GATK3 ------------------
echo ""
echo "[$(timestamp)] Installing GATK3 ..."
gatk3_dir="$build_dir/GATK3"
if [ -z $(check_installed $gatk3_dir) ]; then
    cd $build_dir
    clean "$build_dir/GATK3"
    echo "Create GATK3 folder for users' manual installation"
    mkdir GATK3
    cd GATK3
    # download "https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/${SRA_VERSION}/GenomeAnalysisTK.jar" GenomeAnalysisTK.jar
    git clone $GATK3_DOWNLOAD_URL
    mv ./GATK3_Archive/gatk3.jar.gz .
    gunzip gatk3.jar.gz
    mv gatk3.jar GenomeAnalysisTK.jar
    chmod 755 GenomeAnalysisTK.jar
    ln -s GenomeAnalysisTK.jar gatk3.jar     
    note_installed $gatk3_dir
fi

# --------------- GEM-Tools -----------------
echo ""
echo "[$(timestamp)] Installing GEM-Tools ..."                                                                                        
gemtools_dir="$build_dir/gemtools-${GEMTOOLS_VERSION}-i3/bin"
if [ -z $(check_installed $gemtools_dir) ]; then
    cd $build_dir
    clean "$build_dir/gemtools-${GEMTOOLS_VERSION}-i3"
    echo "Download GEMTOOLS-v${GEMTOOLS_VERSION}"
    download $GEMTOOLS_DOWNLOAD_URL "gemtools-${GEMTOOLS_VERSION}.tar.gz"
    tar xvzf "gemtools-${GEMTOOLS_VERSION}.tar.gz"
    rm gemtools-${GEMTOOLS_VERSION}.tar.gz
    note_installed $gemtools_dir
fi

# --------------- GATK4 ------------------                                                                                                
echo ""
echo "[$(timestamp)] Installing GATK4 ..."
gatk4_dir="$build_dir/GATK4"
if [ -z $(check_installed $gatk4_dir) ]; then
    cd $build_dir
    clean "$build_dir/GATK4"
    echo "Download GATK4-v${GATK4_VERSION}"
    download $GATK4_DOWNLOAD_URL "gatk-${GATK4_VERSION}.zip"
    unzip gatk-${GATK4_VERSION}.zip
    mv gatk-${GATK4_VERSION} GATK4
    cd GATK4
    ln -s gatk-package-${GATK4_VERSION}-local.jar gatk4.jar
    cd $build_dir
    rm gatk-${GATK4_VERSION}.zip
    note_installed $gatk4_dir
fi

# --------------- Freebayes -----------------                                                                                       
echo ""
echo "[$(timestamp)] Installing Freebayes ..."
freebayes_dir="$build_dir/freebayes_conda_env/bin"
if [ -z $(check_installed $freebayes_dir) ]; then
    cd $build_dir
    clean "$build_dir/freebayes_conda_env"
    echo "Download Freebayes-v${FREEBAYES_VERSION}"
    $miniconda3_dir/conda create -y -p $build_dir/freebayes_conda_env python=${PYTHON3_VERSION}
    source $miniconda3_dir/activate $build_dir/freebayes_conda_env
    if [[ "$mainland_china_installation" == "yes" ]]; then
        $miniconda3_dir/conda install -y freebayes=${FREEBAYES_VERSION}
    else
        $miniconda3_dir/conda install -y -c bioconda freebayes=${FREEBAYES_VERSION}
    fi
    source $miniconda3_dir/deactivate
    note_installed $freebayes_dir
fi


# # --------------- Freebayes -----------------
# echo ""
# echo "[$(timestamp)] Installing Freebayes ..."
# freebayes_dir="$build_dir/freebayes-${FREEBAYES_VERSION}/bin"
# if [ -z $(check_installed $freebayes_dir) ]; then
#     cd $build_dir
#     clean "$build_dir/freebayes-${FREEBAYES_VERSION}"
#     echo "Download Freebayes-v${FREEBAYES_VERSION}"
#     download $FREEBAYES_SOURCE_DOWNLOAD_URL "freebayes-${FREEBAYES_VERSION}.tar.gz"
#     tar -xzf freebayes-${FREEBAYES_VERSION}.tar.gz
#     cd freebayes-${FREEBAYES_VERSION}
#     #cd bin
#     #chmod u+x freebayes-${FREEBAYES_VERSION}-linux-static-AMD64
#     #ln -s freebayes-${FREEBAYES_VERSION}-linux-static-AMD64 freebayes
#     cd $build
#     if [ -f freebayes-${FREEBAYES_VERSION}.tar.gz ]; then 
# 	rm freebayes-${FREEBAYES_VERSION}.tar.gz
#     fi
#     note_installed $freebayes_dir
# fi

# --------------- Xatlas -----------------                                                                                       
echo ""
echo "[$(timestamp)] Installing Xatlas ..."
xatlas_dir="$build_dir/xatlas_conda_env/bin"
if [ -z $(check_installed $xatlas_dir) ]; then
    cd $build_dir
    clean "$build_dir/xatlas_conda_env"
    echo "Download Xatlas-v${XATLAS_VERSION}"
    $miniconda3_dir/conda create -y -p $build_dir/xatlas_conda_env python=${PYTHON3_VERSION}
    source $miniconda3_dir/activate $build_dir/xatlas_conda_env
    if [[ "$mainland_china_installation" == "yes" ]]; then
        $miniconda3_dir/conda install -y xatlas=${XATLAS_VERSION}
    else
        $miniconda3_dir/conda install -y -c bioconda xatlas=${XATLAS_VERSION}
    fi
    source $miniconda3_dir/deactivate
    note_installed $xatlas_dir
fi


# --------------- Vcflib -----------------                                                                                       
echo ""
echo "[$(timestamp)] Installing Vcflib ..."
vcflib_dir="$build_dir/vcflib_conda_env/bin"
if [ -z $(check_installed $vcflib_dir) ]; then
    cd $build_dir
    clean "$build_dir/vcflib_conda_env"
    echo "Download Vcflib-v${VCFLIB_VERSION}"
    $miniconda3_dir/conda create -y -p $build_dir/vcflib_conda_env python=${PYTHON3_VERSION}
    source $miniconda3_dir/activate $build_dir/vcflib_conda_env
    if [[ "$mainland_china_installation" == "yes" ]]; then
        $miniconda3_dir/conda install -y vcflib=${VCFLIB_VERSION}
    else
        $miniconda3_dir/conda install -y -c bioconda vcflib=${VCFLIB_VERSION}
    fi
    source $miniconda3_dir/deactivate
    note_installed $vcflib_dir
fi


# --------------- DeepVariant -----------------   
#echo ""
#echo "[$(timestamp)] Installing DeepVariant ..."
#deepvariant_dir="$build_dir/deepvariant_conda_env/bin"
#if [ -z $(check_installed $deepvariant_dir) ]; then
#    cd $build_dir
#    clean "$build_dir/deepvariant_conda_env"
#    echo "Download Deepvariant-v${DEEPVARIANT_VERSION}"
#    $miniconda3_dir/conda create -y -p $build_dir/deepvariant_conda_env python=3.6 # python=${PYTHON3_VERSION}
#    source $miniconda3_dir/activate $build_dir/deepvariant_conda_env
#    if [[ "$mainland_china_installation" == "yes" ]]; then
#	$miniconda3_dir/conda install -y deepvariant=${DEEPVARIANT_VERSION}
#    else
#	$miniconda3_dir/conda install -y -c bioconda deepvariant=${DEEPVARIANT_VERSION}
#    fi
#    source $miniconda3_dir/deactivate
#    note_installed $deepvariant_dir
#fi

# --------------- Clair3 -----------------   
echo ""
echo "[$(timestamp)] Installing Clair3 ..."
clair3_dir="$build_dir/clair3_conda_env/bin"
if [ -z $(check_installed $clair3_dir) ]; then
    cd $build_dir
    clean "$build_dir/clair3_conda_env"
    echo "Download Clair3-v${CLAIR3_VERSION}"
    $miniconda3_dir/conda create -y -p $build_dir/clair3_conda_env python=3.9 # python=${PYTHON3_VERSION}
    source $miniconda3_dir/activate $build_dir/clair3_conda_env
    if [[ "$mainland_china_installation" == "yes" ]]; then
	$miniconda3_dir/conda install -y clair3=${CLAIR3_VERSION} clair3-illumina=${CLAIR3ILLUMINA_VERSION}
    else
	$miniconda3_dir/conda install -y -c bioconda clair3=${CLAIR3_VERSION} clair3-illumina=${CLAIR3ILLUMINA_VERSION}
    fi
    cd $build_dir/clair3_conda_env
    source $miniconda3_dir/deactivate
    note_installed $clair3_dir
fi

# --------------- longshot -----------------
echo ""
echo "[$(timestamp)] Installing longshot ..."   
longshot_dir="$build_dir/longshot_conda_env/bin"
if [ -z $(check_installed $longshot_dir) ]; then
    cd $build_dir
    clean "$build_dir/longshot_conda_env"
    echo "Download longshot-v${LONGSHOT_VERSION}"
    $miniconda3_dir/conda create -y -p $build_dir/longshot_conda_env python=${PYTHON3_VERSION}
    source $miniconda3_dir/activate $build_dir/longshot_conda_env
    if [[ "$mainland_china_installation" == "yes" ]]; then
	$miniconda3_dir/conda install -y longshot=${LONGSHOT_VERSION}
    else
	$miniconda3_dir/conda install -y -c bioconda longshot=${LONGSHOT_VERSION}
    fi
    source $miniconda3_dir/deactivate
    note_installed $longshot_dir    
fi


# --------------- whatshap -----------------
echo ""
echo "[$(timestamp)] Installing whatshap ..."   
whatshap_dir="$build_dir/whatshap_conda_env/bin"
if [ -z $(check_installed $whatshap_dir) ]; then
    cd $build_dir
    clean "$build_dir/whatshap_conda_env"
    echo "Download whatshap-v${WHATSHAP_VERSION}"
    $miniconda3_dir/conda create -y -p $build_dir/whatshap_conda_env python=${PYTHON3_VERSION}
    source $miniconda3_dir/activate $build_dir/whatshap_conda_env
    if [[ "$mainland_china_installation" == "yes" ]]; then
	$miniconda3_dir/conda install -y whatshap=${WHATSHAP_VERSION} nomkl
    else
	$miniconda3_dir/conda install -y -c bioconda whatshap=${WHATSHAP_VERSION} nomkl
    fi
    source $miniconda3_dir/deactivate
    note_installed $whatshap_dir    
fi

# --------------- cnvkit -----------------
echo ""
echo "[$(timestamp)] Installing cnvkit ..."   
cnvkit_dir="$build_dir/cnvkit_conda_env/bin"
if [ -z $(check_installed $cnvkit_dir) ]; then
    cd $build_dir
    clean "$build_dir/cnvkit_conda_env"
    echo "Download cnvkit-v${CNVKIT_VERSION}"
    $miniconda3_dir/conda create -y -p $build_dir/cnvkit_conda_env python=${PYTHON3_VERSION}
    source $miniconda3_dir/activate $build_dir/cnvkit_conda_env
    if [[ "$mainland_china_installation" == "yes" ]]; then
	$miniconda3_dir/conda install -y cnvkit=${CNVKIT_VERSION} nomkl
    else
	$miniconda3_dir/conda install -y -c bioconda cnvkit=${CNVKIT_VERSION} nomkl
    fi
    source $miniconda3_dir/deactivate
    note_installed $cnvkit_dir    
fi

# ----------------- FREEC --------------------
echo ""
echo "[$(timestamp)] Installing FREEC ..."
freec_dir="$build_dir/FREEC-${FREEC_VERSION}"
if [ -z $(check_installed $freec_dir) ]; then
    cd $build_dir
    clean "$build_dir/FREEC-${FREEC_VERSION}"
    echo "Download FREEC-v${FREEC_VERSION}"
    download $FREEC_DOWNLOAD_URL "FREEC-${FREEC_VERSION}.tar.gz"
    tar xvzf FREEC-${FREEC_VERSION}.tar.gz
    cd $freec_dir/src
    make -j $MAKE_JOBS
    cp freec $freec_dir
    cd $freec_dir/scripts
    chmod 755 *
    cp * $freec_dir
    cd $build_dir
    rm FREEC-${FREEC_VERSION}.tar.gz
    note_installed $freec_dir
fi

# ------------------ Manta ---------------------
echo ""
echo "[$(timestamp)] Installing Manta ..."
manta_dir="$build_dir/manta-${MANTA_VERSION}.centos6_x86_64/bin"
if [ -z $(check_installed $manta_dir) ]; then
    cd $build_dir
    clean "$build_dir/manta-${MANTA_VERSION}.centos6_x86_64"
    echo "Download Manta-v${MANTA_VERSION}"
    download $MANTA_DOWNLOAD_URL "manta-${MANTA_VERSION}.centos6_x86_64.tar.bz2"
    tar xvjf manta-${MANTA_VERSION}.centos6_x86_64.tar.bz2
    rm manta-${MANTA_VERSION}.centos6_x86_64.tar.bz2
    note_installed $manta_dir
fi

# ------------------ SVABA ---------------------
echo ""
echo "[$(timestamp)] Installing SVABA ..."
svaba_dir="$build_dir/svaba_conda_env/bin"
if [ -z $(check_installed $svaba_dir) ]; then
    cd $build_dir
    clean "$build_dir/svaba_conda_env"
    echo "Download SVABA-v${SVABA_VERSION}"
    $miniconda3_dir/conda create -y -p $build_dir/svaba_conda_env python=${PYTHON3_VERSION} 
    source $miniconda3_dir/activate $build_dir/svaba_conda_env
    if [[ "$mainland_china_installation" == "yes" ]]; then
    	$miniconda3_dir/conda install -y svaba=${SVABA_VERSION}
    else
    	$miniconda3_dir/conda install -y -c bioconda svaba=${SVABA_VERSION}
    fi
    source $miniconda3_dir/deactivate
    note_installed $svaba_dir
fi

# ------------------ Delly ---------------------
echo ""
echo "[$(timestamp)] Installing Delly ..."
delly_dir="$build_dir/delly-${DELLY_VERSION}"
if [ -z $(check_installed $delly_dir) ]; then
    cd $build_dir
    clean "$build_dir/delly-${DELLY_VERSION}"
    echo "Download Delly-v${DELLY_VERSION}"
    mkdir delly-${DELLY_VERSION}
    cd delly-${DELLY_VERSION}
    download $DELLY_DOWNLOAD_URL "delly_v${DELLY_VERSION}_parallel_linux_x86_64bit"
    chmod 755 delly_v${DELLY_VERSION}_parallel_linux_x86_64bit
    ln -s delly_v${DELLY_VERSION}_parallel_linux_x86_64bit delly
    note_installed $delly_dir
fi

# ------------------ SVIM ---------------------
echo ""
echo "[$(timestamp)] Installing SVIM ..."
svim_dir="$build_dir/svim_conda_env/bin"
if [ -z $(check_installed $svim_dir) ]; then
    cd $build_dir
    clean "$build_dir/svim_conda_env"
    echo "Download SVIM-v${SVIM_VERSION}"
    $miniconda3_dir/conda create -y -p $build_dir/svim_conda_env python=${PYTHON3_VERSION} 
    source $miniconda3_dir/activate $build_dir/svim_conda_env
    # if [[ "$mainland_china_installation" == "yes" ]]; then
    # 	$miniconda3_dir/conda install -y svim=${SVIM_VERSION}
    # else
    # 	$miniconda3_dir/conda install -y -c bioconda svim=${SVIM_VERSION}
    # fi
    $miniconda3_dir/pip3 install svim==${SVIM_VERSION}
    source $miniconda3_dir/deactivate
    note_installed $svim_dir
fi

# ------------------ SNIFFLES ---------------------
echo ""
echo "[$(timestamp)] Installing SNIFFLES ..."
sniffles_dir="$build_dir/sniffles_conda_env/bin"
if [ -z $(check_installed $sniffles_dir) ]; then
    cd $build_dir
    clean "$build_dir/sniffles_conda_env"
    echo "Download SNIFFLES-v${SNIFFLES_VERSION}"
    $miniconda3_dir/conda create -y -p $build_dir/sniffles_conda_env python=${PYTHON3_VERSION} 
    source $miniconda3_dir/activate $build_dir/sniffles_conda_env
    if [[ "$mainland_china_installation" == "yes" ]]; then
    	$miniconda3_dir/conda install -y sniffles=${SNIFFLES_VERSION}
    else
    	$miniconda3_dir/conda install -y -c bioconda sniffles=${SNIFFLES_VERSION}
    fi
    source $miniconda3_dir/deactivate
    note_installed $sniffles_dir
fi

# # ------------------ Sniffles ---------------------
# echo ""
# echo "[$(timestamp)] Installing Sniffles ..."
# #sniffles_dir="$build_dir/Sniffles-${SNIFFLES_VERSION}/bin/sniffles-core-${SNIFFLES_VERSION}"
# sniffles_dir="$build_dir/Sniffles/bin"
# if [ -z $(check_installed $sniffles_dir) ]; then
#     cd $build_dir
#     clean "$build_dir/Sniffles-${SNIFFLES_VERSION}"
#     echo "Download Sniffles-v${SNIFFLES_VERSION}"
#     git clone $SNIFFLES_DOWNLOAD_URL 
#     cd Sniffles
#     git checkout -f -q $SNIFFLES_GITHUB_COMMIT_VERSION
#     # download $SNIFFLES_DOWNLOAD_URL "Sniffles-${SNIFFLES_VERSION}.tar.gz"
#     # tar xvzf Sniffles-${SNIFFLES_VERSION}.tar.gz
#     # cd $build_dir/Sniffles-${SNIFFLES_VERSION}
#     mkdir -p build
#     cd build
#     cmake ..
#     make -j $MAKE_JOBS
#     cd $build_dir
#     rm Sniffles-${SNIFFLES_VERSION}.tar.gz
#     note_installed $sniffles_dir
# fi

# ------------------ PBSV ---------------------
echo ""
echo "[$(timestamp)] Installing PBSV ..."
pbsv_dir="$build_dir/pbsv_conda_env/bin"
if [ -z $(check_installed $pbsv_dir) ]; then
    cd $build_dir
    clean "$build_dir/pbsv_conda_env"
    $miniconda3_dir/conda create -y -p $build_dir/pbsv_conda_env python=${PYTHON3_VERSION}
    source $miniconda3_dir/activate $build_dir/pbsv_conda_env
    if [[ "$mainland_china_installation" == "yes" ]]; then
	$miniconda3_dir/conda install -y pbsv=${PBSV_VERSION}
    else
	$miniconda3_dir/conda install -y -c bioconda pbsv=${PBSV_VERSION}
    fi
    source $miniconda3_dir/deactivate
    note_installed $pbsv_dir
fi

# ------------------ PICKY ---------------------
echo ""
echo "[$(timestamp)] Installing Picky ..."   
picky_dir="$build_dir/Picky/src"
if [ -z $(check_installed $picky_dir) ]; then
    cd $build_dir
    clean "$build_dir/Picky"
    echo "Download PICKY-v${PICKY_VERSION}"
    git clone $PICKY_DOWNLOAD_URL
    cd $build_dir/Picky
    git checkout -f -q $PICKY_GITHUB_COMMIT_VERSION
    note_installed $picky_dir
fi

# ------------------ NANOSV ---------------------
echo ""
echo "[$(timestamp)] Installing NanoSV ..."
nanosv_dir="$build_dir/nanosv_conda_env/bin"
if [ -z $(check_installed $nanosv_dir) ]; then
    cd $build_dir
    clean "$build_dir/nanosv_conda_env"
    $miniconda3_dir/conda create -y -p $build_dir/nanosv_conda_env python=${PYTHON3_VERSION}
    source $miniconda3_dir/activate $build_dir/nanosv_conda_env
    if [[ "$mainland_china_installation" == "yes" ]]; then
	$miniconda3_dir/conda install -y nanosv=${NANOSV_VERSION}
    else
	$miniconda3_dir/conda install -y -c bioconda nanosv=${NANOSV_VERSION}
    fi
    source $miniconda3_dir/deactivate
    note_installed $nanosv_dir
fi

# ------------------ NANOVAR ---------------------
echo ""
echo "[$(timestamp)] Installing NanoVar ..."
nanovar_dir="$build_dir/nanovar_conda_env/bin"
if [ -z $(check_installed $nanovar_dir) ]; then
    cd $build_dir
    clean "$build_dir/nanovar_conda_env"
    $miniconda3_dir/conda create -y -p $build_dir/nanovar_conda_env python=${PYTHON3_VERSION}
    source $miniconda3_dir/activate $build_dir/nanovar_conda_env
    $miniconda3_dir/conda config --set channel_priority disabled
    if [[ "$mainland_china_installation" == "yes" ]]; then
	$miniconda3_dir/conda install -y nanovar=${NANOVAR_VERSION}
    else
	$miniconda3_dir/conda install -y -c bioconda nanovar=${NANOVAR_VERSION}
    fi
    $miniconda3_dir/conda config --set channel_priority flexible
    source $miniconda3_dir/deactivate
    note_installed $nanovar_dir
fi

# ------------------ cuteSV ---------------------
echo ""
echo "[$(timestamp)] Installing cuteSV ..."
cutesv_dir="$build_dir/cutesv_conda_env/bin"
if [ -z $(check_installed $cutesv_dir) ]; then
    cd $build_dir
    clean "$build_dir/cutesv_conda_env"
    $miniconda3_dir/conda create -y -p $build_dir/cutesv_conda_env python=${PYTHON3_VERSION}
    source $miniconda3_dir/activate $build_dir/cutesv_conda_env
    if [[ "$mainland_china_installation" == "yes" ]]; then
	$miniconda3_dir/conda install -y cutesv=${CUTESV_VERSION}
    else
	$miniconda3_dir/conda install -y -c bioconda cutesv=${CUTESV_VERSION}
    fi
    source $miniconda3_dir/deactivate
    note_installed $cutesv_dir
fi


# ------------------ DeBreak ---------------------
echo ""
echo "[$(timestamp)] Installing debreak ..."
debreak_dir="$build_dir/debreak_conda_env/bin"
if [ -z $(check_installed $debreak_dir) ]; then
    cd $build_dir
    clean "$build_dir/debreak_conda_env"
    $miniconda3_dir/conda create -y -p $build_dir/debreak_conda_env python=${PYTHON3_VERSION}
    source $miniconda3_dir/activate $build_dir/debreak_conda_env
    if [[ "$mainland_china_installation" == "yes" ]]; then
	$miniconda3_dir/conda install -y debreak=${DEBREAK_VERSION}
    else
	$miniconda3_dir/conda install -y -c bioconda debreak=${DEBREAK_VERSION}
    fi
    source $miniconda3_dir/deactivate
    note_installed $debreak_dir
fi

# ------------------ JasmineSV ---------------------
echo ""
echo "[$(timestamp)] Installing JasmineSV ..."
jasminesv_dir="$build_dir/jasminesv_conda_env/bin"
if [ -z $(check_installed $jasminesv_dir) ]; then
    cd $build_dir
    clean "$build_dir/jasminesv_conda_env"
    $miniconda3_dir/conda create -y -p $build_dir/jasminesv_conda_env python=${PYTHON3_VERSION}
    source $miniconda3_dir/activate $build_dir/jasminesv_conda_env
    if [[ "$mainland_china_installation" == "yes" ]]; then
	$miniconda3_dir/conda install -y jasminesv=${JASMINESV_VERSION}
    else
	$miniconda3_dir/conda install -y -c bioconda jasminesv=${JASMINESV_VERSION}
    fi
    source $miniconda3_dir/deactivate
    note_installed $jasminesv_dir
fi

# ------------------ Usable_VCF ---------------------
echo ""
echo "[$(timestamp)] Installing Usablevcf ..."
usablevcf_dir="$build_dir/usablevcf_conda_env/bin"
if [ -z $(check_installed $usablevcf_dir) ]; then
    cd $build_dir
    clean "$build_dir/usablevcf_conda_env"
    $miniconda3_dir/conda create -y -p $build_dir/usablevcf_conda_env python=${PYTHON3_VERSION}
    source $miniconda3_dir/activate $build_dir/usablevcf_conda_env
    if [[ "$mainland_china_installation" == "yes" ]]; then
	$miniconda3_dir/conda install -y bcftools htslib vcftools pysam pyvcf 
    else
	$miniconda3_dir/conda install -y -c bioconda bcftools htslib vcftools pysam pyvcf 
    fi
    cd $usablevcf_dir
    git clone https://github.com/ACEnglish/usable_vcf.git
    cp ./usable_vcf/usable_vcf.py .
    source $miniconda3_dir/deactivate
    note_installed $usablevcf_dir
fi

# ------------------ Duphold ---------------------
echo ""
echo "[$(timestamp)] Installing duphold ..."
duphold_dir="$build_dir/Duphold-${DUPHOLD_VERSION}"
if [ -z $(check_installed $duphold_dir) ]; then
    cd $build_dir
    clean "$build_dir/Duphold-${DUPHOLD_VERSION}"
    echo "Download Duphold-v${DUPHOLD_VERSION}"
    mkdir Duphold-${DUPHOLD_VERSION}
    cd Duphold-${DUPHOLD_VERSION}
    download $DUPHOLD_DOWNLOAD_URL "duphold"
    chmod u+x duphold
    note_installed $duphold_dir
fi

# ------------------ Truvari ---------------------
echo ""
echo "[$(timestamp)] Installing Truvari ..."
truvari_dir="$build_dir/truvari_conda_env/bin"
if [ -z $(check_installed $truvari_dir) ]; then
    cd $build_dir
    clean "$build_dir/truvari_conda_env"
    $miniconda3_dir/conda create -y -p $build_dir/truvari_conda_env python=${PYTHON3_VERSION}
    source $miniconda3_dir/activate $build_dir/truvari_conda_env
    $miniconda3_dir/python -m pip install truvari==${TRUVARI_VERSION}
    source $miniconda3_dir/deactivate
    note_installed $truvari_dir
fi

# ------------------ Sansa ---------------------
echo ""
echo "[$(timestamp)] Installing sansa ..."
sansa_dir="$build_dir/Sansa-${SANSA_VERSION}"
if [ -z $(check_installed $sansa_dir) ]; then
    cd $build_dir
    clean "$build_dir/Sansa-${SANSA_VERSION}"
    echo "Download Sansa-v${SANSA_VERSION}"
    mkdir Sansa-${SANSA_VERSION}
    cd Sansa-${SANSA_VERSION}
    download $SANSA_DOWNLOAD_URL sansa
    chmod u+x sansa
    note_installed $sansa_dir
fi

# ------------------ Samplot ---------------------
echo ""
echo "[$(timestamp)] Installing Samplot ..."
samplot_dir="$build_dir/samplot_conda_env/bin"
if [ -z $(check_installed $samplot_dir) ]; then
    cd $build_dir
    clean "$build_dir/samplot_conda_env"
    $miniconda3_dir/conda create -y -p $build_dir/samplot_conda_env python=${PYTHON3_VERSION}
    source $miniconda3_dir/activate $build_dir/samplot_conda_env
    if [[ "$mainland_china_installation" == "yes" ]]; then
	$miniconda3_dir/conda install -y samplot=${SAMPLOT_VERSION}
    else
	$miniconda3_dir/conda install -y -c bioconda samplot=${SAMPLOT_VERSION}
    fi
    source $miniconda3_dir/deactivate
    note_installed $samplot_dir
fi

# ------------------ VT ---------------------
echo ""
echo "[$(timestamp)] Installing VT ..."
vt_dir="$build_dir/vt"
if [ -z $(check_installed $vt_dir) ]; then
    cd $build_dir
    clean "$build_dir/vt"
    echo "Download VT-v${VT_VERSION}"
    git clone $VT_DOWNLOAD_URL
    cd $vt_dir
    git checkout -f -q $VT_GITHUB_COMMIT_VERSION
    make -j $MAKE_JOBS
    make test
    note_installed $vt_dir
fi

# --------------- bcftools ------------------
echo ""
echo "[$(timestamp)] Installing bcftools ..."
bcftools_dir="$build_dir/bcftools-${BCFTOOLS_VERSION}"
if [ -z $(check_installed $bcftools_dir) ]; then
    cd $build_dir
    clean "$build_dir/bcftools-${BCFTOOLS_VERSION}"
    echo "Download bcftools-v${BCFTOOLS_VERSION}"
    download $BCFTOOLS_DOWNLOAD_URL "bcftools-${BCFTOOLS_VERSION}.tar.bz2"
    tar xvjf bcftools-${BCFTOOLS_VERSION}.tar.bz2
    cd $bcftools_dir
    ./configure
    make -j $MAKE_JOBS
    cd $build_dir
    rm bcftools-${BCFTOOLS_VERSION}.tar.bz2
    note_installed $bcftools_dir
fi

# --------------- bedtools ------------------
echo ""
echo "[$(timestamp)] Installing bedtools ..."
bedtools_dir="$build_dir/bedtools2/bin"
if [ -z $(check_installed $bedtools_dir) ]; then
    cd $build_dir
    clean "$build_dir/bedtools2"
    echo "Download bedtools-v${BEDTOOLS_VERSION}"
    download $BEDTOOLS_DOWNLOAD_URL "bedtools-${BEDTOOLS_VERSION}.tar.gz"
    tar xvzf bedtools-${BEDTOOLS_VERSION}.tar.gz
    cd "$build_dir/bedtools2"
    make -j $MAKE_JOBS
    cd $build_dir
    rm bedtools-${BEDTOOLS_VERSION}.tar.gz
    note_installed $bedtools_dir
fi

# --------------- ncbi-blast+ ------------------
echo ""
echo "[$(timestamp)] Installing ncbi-blast+ ..."
blast_dir="$build_dir/ncbi-blast-${BLAST_VERSION}+/bin"
windowmasker_dir="$build_dir/ncbi-blast-${BLAST_VERSION}+/bin"
if [ -z $(check_installed $blast_dir) ]; then
    cd $build_dir
    clean "$build_dir/ncbi-blast-${BLAST_VERSION}+"
    echo "Download ncbi-blast-v${BLAST_VERSION}"
    download $BLAST_DOWNLOAD_URL "ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz"
    tar xvzf ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz
    rm ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz
    note_installed $blast_dir
fi

# --------------- ncbi-rmblast ------------------
echo ""
echo "[$(timestamp)] Installing ncbi-rmblast ..."
rmblast_dir="$build_dir/ncbi-rmblastn-${RMBLAST_VERSION}/bin"
if [ -z $(check_installed $rmblast_dir) ]; then
    cd $build_dir
    clean "$build_dir/ncbi-rmblastn-${RMBLAST_VERSION}"
    echo "Download ncbi-rmblastn-v${BLAST_VERSION}"
    download $RMBLAST_DOWNLOAD_URL "ncbi-rmblastn-${RMBLAST_VERSION}-x64-linux.tar.gz"
    tar xvzf ncbi-rmblastn-${RMBLAST_VERSION}-x64-linux.tar.gz
    # copy rmblastn binary file to ncbi-blast+ directory for easy RepeatMasker configuration
    cp $rmblast_dir/rmblastn $blast_dir
    rm ncbi-rmblastn-${RMBLAST_VERSION}-x64-linux.tar.gz
    note_installed $rmblast_dir
fi

# --------------- ensembl-vep -----------------
echo ""
echo "[$(timestamp)] Installing ensembl-vep ..."   
vep_dir="$build_dir/ensembl-vep_conda_env/bin"
if [ -z $(check_installed $vep_dir) ]; then
    cd $build_dir
    clean "$build_dir/ensembl-vep_conda_env"
    echo "Download ensembl-vep-v${VEP_VERSION}"
    $miniconda3_dir/conda create -y -p $build_dir/ensembl-vep_conda_env python=${PYTHON3_VERSION}
    source $miniconda3_dir/activate $build_dir/ensembl-vep_conda_env
    if [[ "$mainland_china_installation" == "yes" ]]; then
        $miniconda3_dir/conda install -y ensembl-vep=${VEP_VERSION}
    else
        $miniconda3_dir/conda install -y -c bioconda ensembl-vep=${VEP_VERSION}
        $miniconda3_dir/conda install -y -c conda-forge perl-compress-raw-zlib=2.202
    fi
    source $miniconda3_dir/deactivate
    note_installed $vep_dir    
fi

# --------------- parallel ------------------
echo ""
echo "[$(timestamp)] Installing parallel ..."
parallel_dir="$build_dir/parallel-${PARALLEL_VERSION}/bin"
if [ -z $(check_installed $parallel_dir) ]; then
    cd $build_dir
    clean "$build_dir/parallel-${PARALLEL_VERSION}"
    echo "Download parallel"
    download $PARALLEL_DOWNLOAD_URL "parallel_v${PARALLEL_VERSION}.tar.bz2"
    tar xvjf parallel_v${PARALLEL_VERSION}.tar.bz2
    cd parallel-${PARALLEL_VERSION}
    ./configure --prefix="$build_dir/parallel-${PARALLEL_VERSION}"
    make -j $MAKE_JOBS
    make install
    parallel_dir="$build_dir/parallel-${PARALLEL_VERSION}/bin"
    cd ..
    rm parallel_v${PARALLEL_VERSION}.tar.bz2
    note_installed $parallel_dir
fi

$miniconda3_dir/conda config --set ssl_verify yes

# Configure executable paths

cd $VARATHON_HOME
echo ""
echo "Configuring executable paths ..."

echo "export VARATHON_HOME=${VARATHON_HOME}" > env.sh
echo "export build_dir=${build_dir}" >> env.sh
echo "export PERL5LIB=${PERL5LIB}" >> env.sh 
echo "export R_LIBS=${R_LIBS}" >> env.sh
echo "export cpanm_dir=${cpanm_dir}" >> env.sh
#echo "export miniconda2_dir=${miniconda2_dir}" >> env.sh
echo "export miniconda3_dir=${miniconda3_dir}" >> env.sh
echo "export simuG_dir=${simuG_dir}" >> env.sh
echo "export nanoplot_dir=${nanoplot_dir}" >> env.sh
echo "export sra_dir=${sra_dir}" >> env.sh
echo "export art_dir=${art_dir}" >> env.sh
echo "export pbsim3_dir=${pbsim3_dir}" >> env.sh
#echo "export simlord_dir=${simlord_dir}" >> env.sh
#echo "export nanosim_dir=${nanosim_dir}" >> env.sh
#echo "export deepsimulator_dir=${deepsimulator_dir}" >> env.sh
echo "export trimmomatic_dir=${trimmomatic_dir}" >> env.sh
echo "export porechop_dir=${porechop_dir}" >> env.sh
echo "export filtlong_dir=${filtlong_dir}" >> env.sh
echo "export bwa_dir=${bwa_dir}" >> env.sh
echo "export bwamem2_dir=${bwamem2_dir}" >> env.sh
echo "export last_dir=${last_dir}" >> env.sh
echo "export ngmlr_dir=${ngmlr_dir}" >> env.sh
echo "export minimap2_dir=${minimap2_dir}" >> env.sh
echo "export winnowmap_dir=${winnowmap_dir}" >> env.sh
echo "export ira_dir=${ira_dir}" >> env.sh
echo "export pbmm2_dir=${pbmm2_dir}" >> env.sh
echo "export graphmap_dir=${graphmap_dir}" >> env.sh
echo "export graphmap2_dir=${graphmap2_dir}" >> env.sh
echo "export samtools_dir=${samtools_dir}" >> env.sh
echo "export htslib_dir=${htslib_dir}" >> env.sh
echo "export tabix_dir=${tabix_dir}" >> env.sh
echo "export picard_dir=${picard_dir}" >> env.sh
echo "export gatk3_dir=${gatk3_dir}" >> env.sh
echo "export gemtools_dir=${gemtools_dir}" >> env.sh
echo "export gatk4_dir=${gatk4_dir}" >> env.sh
echo "export xatlas_dir=${xatlas_dir}" >> env.sh
echo "export deepvariant_dir=${deepvariant_dir}" >> env.sh
echo "export freebayes_dir=${freebayes_dir}" >> env.sh
echo "export vcflib_dir=${vcflib_dir}" >> env.sh
echo "export vt_dir=${vt_dir}" >> env.sh
echo "export clair3_dir=${clair3_dir}" >> env.sh
echo "export longshot_dir=${longshot_dir}" >> env.sh
echo "export whatshap_dir=${whatshap_dir}" >> env.sh
echo "export cnvkit_dir=${cnvkit_dir}" >> env.sh
echo "export freec_dir=${freec_dir}" >> env.sh
echo "export manta_dir=${manta_dir}" >> env.sh
echo "export svaba_dir=${svaba_dir}" >> env.sh
echo "export delly_dir=${delly_dir}" >> env.sh
echo "export svim_dir=${svim_dir}" >> env.sh
echo "export sniffles_dir=${sniffles_dir}" >> env.sh
echo "export pbsv_dir=${pbsv_dir}" >> env.sh
echo "export picky_dir=${picky_dir}" >> env.sh
echo "export nanosv_dir=${nanosv_dir}" >> env.sh
echo "export nanovar_dir=${nanovar_dir}" >> env.sh
echo "export cutesv_dir=${cutesv_dir}" >> env.sh
echo "export debreak_dir=${debreak_dir}" >> env.sh
echo "export duphold_dir=${duphold_dir}" >> env.sh
echo "export jasminesv_dir=${jasminesv_dir}" >> env.sh
echo "export usablevcf_dir=${usablevcf_dir}" >> env.sh
echo "export truvari_dir=${truvari_dir}" >> env.sh
echo "export samplot_dir=${samplot_dir}" >> env.sh
echo "export bcftools_dir=${bcftools_dir}" >> env.sh
echo "export vep_dir=${vep_dir}" >> env.sh
echo "export parallel_dir=${parallel_dir}" >> env.sh
echo "export blast_dir=${blast_dir}" >> env.sh
echo "export rmblast_dir=${rmblast_dir}" >> env.sh
echo "export windowmasker_dir=${windowmasker_dir}" >> env.sh
echo "export bedtools_dir=${bedtools_dir}" >> env.sh


# test java configuration: requireds java 1.8 
echo ""
echo "##########################################"
echo "Testing java configuration ..."
echo ""
java_bin=""
if type -p java
then 
    java_bin=$(which java)
    echo "found java executable in PATH: $java_bin"
elif [[ -n "$JAVA_HOME" ]] && [[ -x "$JAVA_HOME/bin/java" ]]
then 
    java_bin="$JAVA_HOME/bin/java"
    echo "found java executable in JAVA_HOME: $java_bin" 
else 
    echo "";
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
    echo "Failed to detect Java installation in the system!"
    echo "Please install java 1.8, which is a dependency of Varathon!\n";
    echo "After the java installation, please manually set the directory path to java 1.8 executable on the last line of the env.sh file generated by this installation script!"
    echo "export java_dir=" >> env.sh
fi  

if [[ -n "$java_bin" ]]
then
    java_version=$("$java_bin" -version 2>&1 | awk -F '"' '/version/ {print $2}')
    echo "detected java_version: $java_version"
    if [ $(tidy_version "$java_version") -eq $(tidy_version "1.8") ]
    then
	java_dir=$(dirname $java_bin)
	echo "export java_dir=${java_dir}" >> env.sh
        echo "You have the correct java version for Varathon! Varathon will take care of the configuration."
    else
	echo "";
	echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
	echo "Your java version is not the version required by Varathon (java v1.8)!"
        echo "Please manually set the directory path to java 1.8 executable on the last line of the env.sh file generated by this installation script!"
	echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!";
	echo "export java_dir=" >> env.sh
    fi
fi

echo ""
# echo "uncompress large supporting files ..."
# gunzip $VARATHON_HOME/data/*.gz

echo "The installation of all dependencies is complete!"
echo ""
echo ""

############################
# checking Bash exit status

if [ $? -eq 0 ]
then
    echo "#######################################################################"
    echo ""
    echo "Varathon message: This bash script has been successfully processed! :)"
    echo ""
    echo "#######################################################################"
    echo ""
    exit 0
fi
############################
