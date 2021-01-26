#!/bin/bash
# last update: 2021.01.20
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
	
download () {
  url=$1
  download_location=$2
  echo "Downloading $url to $download_location"
  wget -c --no-check-certificate $url -O $download_location
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
    echo ""
    echo "Installing LRSDAY build dependencies for Debian/Ubuntu."
    echo "sudo privileges are required and you will be prompted to enter your password"
    sudo apt-get update
    xargs -a debiandeps sudo apt-get install -y
fi

echo ""

while getopts ":hc" opt
do
    case "${opt}" in
        h)
	    echo "Usage:"
	    echo "bash install_dependencies.sh"
            echo "When installing within mainland China, please run this script with the '-c' option >"
	    echo "bash install_dependencies.sh -c";;
        c)
	    echo "Detected the '-c' option >"
	    echo "Set installation location as 'mainland_china'" 
            mainland_china_installation="yes";;
    esac
done
echo "";

# genome simulation
SIMUG_VERSION="1.0.0" # 
SIMUG_GITHUB_COMMIT_VERSION="940b961" # commited on 2020.06.22
SIMUG_DOWNLOAD_URL="https://github.com/yjx1217/simuG.git"

# reads retrieving/simulation/processing
SRA_VERSION="2.9.6" # released on 2019.03.18
SRA_DOWNLOAD_URL="https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/${SRA_VERSION}/sratoolkit.${SRA_VERSION}-centos_linux64.tar.gz"

ART_VERSION="mountrainier2016.06.05" # released on 2016.06.05
ART_DOWNLOAD_URL="https://www.niehs.nih.gov/research/resources/assets/docs/artbin${ART_VERSION}linux64.tgz"

SIMLORD_VERSION="1.0.4" # released on 2018.06.30

DEEPSIMULATOR_VERSION="1.5.0" # released on 2019.06.16 
DEEPSIMULATOR_GITHUB_COMMIT_VERSION="3c867c2" # committed on 2019.11.02
DEEPSIMULATOR_DOWNLOAD_URL="https://github.com/lykaust15/DeepSimulator"

TRIMMOMATIC_VERSION="0.38" # released on 
TRIMMOMATIC_DOWNLOAD_URL="http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-${TRIMMOMATIC_VERSION}.zip"

PORECHOP_VERSION="0.2.4" # 
PORECHOP_GITHUB_COMMIT_VERSION="109e437" # committed on 2018.10.19
PORECHOP_DOWNLOAD_URL="https://github.com/rrwick/Porechop"

FILTLONG_VERSION="0.2.0" #
FILTLONG_GITHUB_COMMIT_VERSION="d1bb46d" # committed on 2018.05.11
FILTLONG_DOWNLOAD_URL="https://github.com/rrwick/Filtlong"

# alignment building/processing
BWA_VERSION="0.7.17" # released on 2017.10.23
BWA_DOWNLOAD_URL="https://github.com/lh3/bwa/releases/download/v${BWA_VERSION}/bwa-${BWA_VERSION}.tar.bz2"

LAST_VERSION="979" # released on 
LAST_DOWNLOAD_URL="http://last.cbrc.jp/last-${LAST_VERSION}.zip"

NGMLR_VERSION="0.2.7" # released on 2018.06.25
NGMLR_DOWNLOAD_URL="https://github.com/philres/ngmlr/releases/download/v${NGMLR_VERSION}/ngmlr-${NGMLR_VERSION}-linux-x86_64.tar.gz"

MINIMAP2_VERSION="2.16" # released on 2019.02.28
MINIMAP2_DOWNLOAD_URL="https://github.com/lh3/minimap2/releases/download/v${MINIMAP2_VERSION}/minimap2-${MINIMAP2_VERSION}_x64-linux.tar.bz2"

PBMM2_VERSION="1.0.0" # released on 2019.03.11
# distributed via conda

GRAPHMAP_VERSION="0.5.2"
GRAPHMAP_GITHUB_COMMIT_VERSION="eb8c75d"
GRAPHMAP_DOWNLOAD_URL="https://github.com/isovic/graphmap"

GRAPHMAP2_VERSION="0.6.4"
GRAPHMAP2_DOWNLOAD_URL="https://github.com/lbcb-sci/graphmap2/releases/download/v0.6.4/prebuild.zip"

SAMTOOLS_VERSION="1.9" # released on 2018.07.18
SAMTOOLS_DOWNLOAD_URL="https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2"

PICARD_VERSION="2.19.0" # released on 2019.03.22
PICARD_DOWNLOAD_URL="https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard.jar"

GATK3_VERSION="3.6-6-g965413b" # 
GATK3_DOWLOAD_URL="https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/${SRA_VERSION}/GenomeAnalysisTK.jar"

GEMTOOLS_VERSION="1.7.1"
GEMTOOLS_DOWNLOAD_URL="http://barnaserver.com/gemtools/releases/GEMTools-static-i3-${GEMTOOLS_VERSION}.tar.gz"

# variant calling
GATK4_VERSION="4.0.11.0" # released on 2018.10.23
GATK4_DOWNLOAD_URL="https://github.com/broadinstitute/gatk/releases/download/${GATK4_VERSION}/gatk-${GATK4_VERSION}.zip"

FREEBAYES_VERSION="1.2.0" # 
FREEBAYES_GITHUB_COMMIT_VERSION="d15209e" # committed on 2019.02.14
FREEBAYES_DOWNLOAD_URL="https://github.com/ekg/freebayes"

CLAIR_VERSION="2.1.0" # released on 2020.04.15

LONGSHOT_VERSION="0.3.4" # released on 2019.05.02

FREEC_VERSION="11.4" # released on 2018.04.28
FREEC_DOWNLOAD_URL="https://github.com/BoevaLab/FREEC/archive/v${FREEC_VERSION}.tar.gz"

MANTA_VERSION="1.5.0" # released on 2018.11.12
MANTA_DOWNLOAD_URL="https://github.com/Illumina/manta/releases/download/v${MANTA_VERSION}/manta-${MANTA_VERSION}.centos6_x86_64.tar.bz2"

DELLY_VERSION="0.8.1" # released on 2019.02.04
DELLY_DOWNLOAD_URL="https://github.com/dellytools/delly/releases/download/v${DELLY_VERSION}/delly_v${DELLY_VERSION}_linux_x86_64bit"

SVIM_VERSION="1.0.0" # released on 2019.04.29
SVIM_GITHUB_COMMIT_VERSION="b72e631" # committed on 2019.05.02

SNIFFLES_VERSION="1.0.11" # released on 2019.01.30
SNIFFLES_GITHUB_COMMIT_VERSION="e5d5150" # committed on 2019.03.30
SNIFFLES_DOWNLOAD_URL="https://github.com/fritzsedlazeck/Sniffles/archive/${SNIFFLES_VERSION}.tar.gz"

PBSV_VERSION="2.2.0" # released on 2019.02.28

PICKY_VERSION="0.2.a" # released on 2018.07.17   
PICKY_GITHUB_COMMIT_VERSION="34b85ac" # committed on 2018.07.17
PICKY_DOWNLOAD_URL="https://github.com/TheJacksonLaboratory/Picky"

NANOSV_VERSION="1.2.3" # released on 2019.04.09
NANOSV_GITHUB_COMMIT_VERSION="c1ae30c" # committed on 2019.04.09
NANOSV_DOWNLOAD_URL="https://github.com/mroosmalen/nanosv"

# variant processing
VT_VERSION="" # we use the github commit version below
VT_GITHUB_COMMIT_VERSION="6686b5c" # committed on 2017.05.09
VT_DOWNLOAD_URL="https://github.com/atks/vt"

HTSLIB_VERSION="1.9" # released on 2018.07.18
HTSLIB_DOWNLOAD_URL="https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2"

BCFTOOLS_VERSION="1.9" # released on 2018.07.18
BCFTOOLS_DOWNLOAD_URL="https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2"

# others
MINICONDA2_VERSION="4.5.11" # released on 2018.09.04
if [[ "$mainland_china_installation" == "no" ]]
then
    MINICONDA2_DOWNLOAD_URL="https://repo.continuum.io/miniconda/Miniconda2-${MINICONDA2_VERSION}-Linux-x86_64.sh"
else
    MINICONDA2_DOWNLOAD_URL="https://mirrors.tuna.tsinghua.edu.cn/anaconda/miniconda/Miniconda2-${MINICONDA2_VERSION}-Linux-x86_64.sh"
fi

MINICONDA3_VERSION="py37_4.8.2" # released on 2020.03.11
if [[ "$mainland_china_installation" == "no" ]]
then
    MINICONDA3_DOWNLOAD_URL="https://repo.anaconda.com/miniconda/Miniconda3-${MINICONDA3_VERSION}-Linux-x86_64.sh"
else
    MINICONDA3_DOWNLOAD_URL="https://mirrors.tuna.tsinghua.edu.cn/anaconda/miniconda/Miniconda3-${MINICONDA3_VERSION}-Linux-x86_64.sh"
fi

BEDTOOLS_VERSION="2.27.1" # released on 2017.12.15
BEDTOOLS_DOWNLOAD_URL="https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLS_VERSION}/bedtools-${BEDTOOLS_VERSION}.tar.gz"

BLAST_VERSION="2.2.31"
BLAST_DOWNLOAD_URL="https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${BLAST_VERSION}/ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz"

RMBLAST_VERSION="2.2.28"
RMBLAST_DOWNLOAD_URL="https://ftp.ncbi.nlm.nih.gov/blast/executables/rmblast/${RMBLAST_VERSION}/ncbi-rmblastn-${RMBLAST_VERSION}-x64-linux.tar.gz"

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

echo ""
echo "[$(timestamp)] Installing Perl modules ..."
cpanm_dir="$build_dir/cpanm"
if [ -z $(check_installed $cpanm_dir) ]; then
    mkdir -p  $cpanm_dir
    cd $cpanm_dir
    wget -c --no-check-certificate -O - https://cpanmin.us/ > cpanm
    chmod +x cpanm
    mkdir perlmods
    $cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed Test::More@1.302086
    $cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed Env@1.04
    $cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed Statistics::Descriptive@3.0612
    $cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed Statistics::Descriptive::Discrete@0.07
    $cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed Math::Random@0.72
    $cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed Math::Round@0.07
    note_installed $cpanm_dir
fi    

echo ""
echo "[$(timestamp)] Installing R libraries ..."
rlib_dir="$build_dir/R_libs"
if [ -z $(check_installed $rlib_dir) ]; then
    mkdir -p $rlib_dir
    cd $rlib_dir
    R_VERSION=$(R --version |head -1 |cut -d " " -f 3)
    R -e "install.packages(\"optparse\", repos=\"http://cran.rstudio.com/\", lib=\"$build_dir/R_libs/\")"
    R -e "install.packages(\"ggplot2\", repos=\"http://cran.rstudio.com/\", lib=\"$build_dir/R_libs/\")"
    R -e "install.packages(\"scales\", repos=\"http://cran.rstudio.com/\", lib=\"$build_dir/R_libs/\")"
    R -e "install.packages(\"viridis\", repos=\"http://cran.rstudio.com/\", lib=\"$build_dir/R_libs/\")"
    if [ $(tidy_version "$R_VERSION") -ge $(tidy_version "3.6.0") ]
    then
	echo "R_VERSION=$R_VERSION, use the new bioconductor installation protocol"
	R -e ".libPaths(\"$build_dir/R_libs/\");install.packages(\"BiocManager\", repos=\"http://cran.rstudio.com/\", lib=\"$build_dir/R_libs/\");BiocManager::install(\"D\
NAcopy\")"
    else
	echo "R_VERSION=$R_VERSION, use the old bioconductor installation protocol"
	R -e ".libPaths(\"$build_dir/R_libs/\");source(\"https://bioconductor.org/biocLite.R\");biocLite(\"DNAcopy\", type = \"source\")"
    fi
    note_installed $rlib_dir
fi

# install dependencies

# ------------- Miniconda2 --------------------
echo ""
echo "[$(timestamp)] Installing miniconda2 ..."
miniconda2_dir="$build_dir/miniconda2/bin"
if [ -z $(check_installed $miniconda2_dir) ]; then
    cd $build_dir
    clean "$build_dir/miniconda2"
    download $MINICONDA2_DOWNLOAD_URL "Miniconda2-${MINICONDA2_VERSION}-Linux-x86_64.sh"
    bash Miniconda2-${MINICONDA2_VERSION}-Linux-x86_64.sh -b -p $build_dir/miniconda2
    if [[ "$mainland_china_installation" == "yes" ]]
    then
	$miniconda2_dir/conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/
	$miniconda2_dir/conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/
	$miniconda2_dir/conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/pytorch/
	$miniconda2_dir/conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/
	$miniconda2_dir/conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
    else 
	$miniconda2_dir/conda config --add channels defaults
	$miniconda2_dir/conda config --add channels bioconda
	$miniconda2_dir/conda config --add channels conda-forge
    fi
    $miniconda2_dir/conda config --set show_channel_urls yes
    $miniconda2_dir/pip install cython==0.29.14 numpy==1.13.1
    rm Miniconda2-${MINICONDA2_VERSION}-Linux-x86_64.sh 
    note_installed $miniconda2_dir
fi

# ------------- Miniconda3 --------------------
echo ""
echo "[$(timestamp)] Installing miniconda3 ..."
miniconda3_dir="$build_dir/miniconda3/bin"
if [ -z $(check_installed $miniconda3_dir) ]; then
    cd $build_dir
    clean "$build_dir/miniconda3"
    download $MINICONDA3_DOWNLOAD_URL "Miniconda3-${MINICONDA3_VERSION}-Linux-x86_64.sh"
    bash Miniconda3-${MINICONDA3_VERSION}-Linux-x86_64.sh -b -p $build_dir/miniconda3
    if [[ "$mainland_china_installation" == "yes" ]]
    then
        $miniconda3_dir/conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/
        $miniconda3_dir/conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/
        $miniconda3_dir/conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/pytorch/
        $miniconda3_dir/conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/
        $miniconda3_dir/conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
    else 
        $miniconda3_dir/conda config --add channels defaults
        $miniconda3_dir/conda config --add channels bioconda
        $miniconda3_dir/conda config --add channels conda-forge
    fi
    $miniconda3_dir/conda config --set show_channel_urls yes
    # $miniconda3_dir/conda install -y -c conda-forge pypy3.7 
    $miniconda3_dir/conda install -y -c conda-forge cython=0.29.21 numpy==1.19.5
    # $miniconda3_dir/pip install cython==0.29.14 numpy==1.13.1
    cd $build_dir
    rm Miniconda3-${MINICONDA3_VERSION}-Linux-x86_64.sh 
    note_installed $miniconda3_dir
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
echo "[$(timestamp)] Installing SRAtoolkit ..."
sra_dir="$build_dir/sratoolkit.${SRA_VERSION}-centos_linux64/bin"
if [ -z $(check_installed $sra_dir) ]; then
    cd $build_dir
    clean "$build_dir/sratoolkit.${SRA_VERSION}-centos_linux64"
    echo "Download SRAtoolkit-v${SRA_VERSION}"
    download $SRA_DOWNLOAD_URL sratoolkit.${SRA_VERSION}-centos_linux64.tar.gz
    tar xvzf sratoolkit.${SRA_VERSION}-centos_linux64.tar.gz
    rm sratoolkit.${SRA_VERSION}-centos_linux64.tar.gz
    note_installed $sra_dir
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

# -------------- SimLoRD ----------------
echo ""
echo "[$(timestamp)] Installing SimLoRD ..."
simlord_dir="$build_dir/simlord_conda_env/bin"
if [ -z $(check_installed $simlord_dir) ]; then
    cd $build_dir
    clean "$build_dir/simlord_conda_env"
    $miniconda3_dir/conda create -y -p $build_dir/simlord_conda_env python=3.7
    source $miniconda3_dir/activate $build_dir/simlord_conda_env
    $miniconda3_dir/conda install -y -c bioconda simlord=${SIMLORD_VERSION} 
    # $miniconda3_dir/pip install pysam==0.15.3
    # $miniconda3_dir/pip install dinopy==2.0.3
    # $miniconda3_dir/pip install "simlord==${SIMLORD_VERSION}"
    source $miniconda3_dir/deactivate
    note_installed $simlord_dir
fi

# -------------- DeepSimulator ----------------
echo ""
echo "[$(timestamp)] Installing DeepSimulator ..."
deepsimulator_dir="$build_dir/DeepSimulator"
if [ -z $(check_installed $deepsimulator_dir) ]; then
    cd $build_dir
    clean "$build_dir/DeepSimulator"
    echo "Download DeepSimulator-v${DEEPSIMULATOR_GITHUB_COMMIT_VERSION}"
    # git clone $DEEPSIMULATOR_DOWNLOAD_URL
    clone $DEEPSIMULATOR_DOWNLOAD_URL
    cd $deepsimulator_dir
    git checkout -f -q $DEEPSIMULATOR_GITHUB_COMMIT_VERSION
    $miniconda2_dir/conda remove --name tensorflow_cdpm --all -y
    $miniconda2_dir/conda create --name tensorflow_cdpm python=2.7 -y
    source $miniconda2_dir/activate tensorflow_cdpm
    $miniconda2_dir/conda install -y -c anaconda scikit-learn=0.20.3
    $miniconda2_dir/pip install numpy==1.13.1
    if [[ "$mainland_china_installation" == "yes" ]]
    then
	$miniconda2_dir/pip install tensorflow==1.2.1 -i https://pypi.tuna.tsinghua.edu.cn/simple/
    else
	$miniconda2_dir/pip install tensorflow==1.2.1
    fi
    $miniconda2_dir/pip install tflearn==0.3.2
    $miniconda2_dir/pip install tqdm==4.19.4
    $miniconda2_dir/pip install scipy==0.18.1
    $miniconda2_dir/pip install h5py==2.7.1
    $miniconda2_dir/pip install biopython==1.74
    source $miniconda2_dir/deactivate
    $miniconda2_dir/conda remove --name basecall --all -y
    $miniconda2_dir/conda create --name basecall python=3.6 -y
    source $miniconda2_dir/activate basecall
    #-> 2. install basecaller
    #--| 2.1 install albacore_2.3.1
    cd base_caller/albacore_2.3.1/
    wget -q https://mirror.oxfordnanoportal.com/software/analysis/ont_albacore-2.3.1-cp36-cp36m-manylinux1_x86_64.whl
    pip install ont_albacore-2.3.1-cp36-cp36m-manylinux1_x86_64.whl
    rm -f ont_albacore-2.3.1-cp36-cp36m-manylinux1_x86_64.whl
    cd ../../
    #--| 2.2 install guppy_3.1.5
    cd base_caller/guppy_3.1.5/
    wget -q https://mirror.oxfordnanoportal.com/software/analysis/ont-guppy-cpu_3.1.5_linux64.tar.gz
    tar xzf ont-guppy-cpu_3.1.5_linux64.tar.gz
    rm -f ont-guppy-cpu_3.1.5_linux64.tar.gz
    cd ../../
    source $miniconda2_dir/deactivate
    cd $deepsimulator_dir
    note_installed $deepsimulator_dir
fi
# --------------- Trimmomatic -----------------
echo ""
echo "[$(timestamp)] Installing Trimmomatic ..."
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
    $miniconda3_dir/conda create -y -p $build_dir/porechop_conda_env python=3.7
    source $miniconda3_dir/activate $build_dir/porechop_conda_env
    $miniconda3_dir/conda install -y -c bioconda porechop=${PORECHOP_VERSION} 
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
echo "[$(timestamp)] Installing BWA ..."
bwa_dir="$build_dir/bwa-${BWA_VERSION}"
if [ -z $(check_installed $bwa_dir) ]; then
    cd $build_dir
    clean "$build_dir/bwa-${BWA_VERSION}"
    echo "Download BWA-v${BWA_VERSION}"
    download $BWA_DOWNLOAD_URL "bwa-${BWA_VERSION}.tar.bz2"
    tar -xjf bwa-${BWA_VERSION}.tar.bz2
    cd $bwa_dir
    make -j $MAKE_JOBS
    cd $build_dir
    rm bwa-${BWA_VERSION}.tar.bz2
    note_installed $bwa_dir
fi

# ------------- LAST -------------------
echo ""
echo "[$(timestamp)] Installing LAST ..."
last_dir="$build_dir/last-${LAST_VERSION}/bin"
if [ -z $(check_installed $last_dir) ]; then
    cd $build_dir
    clean "$build_dir/last-${LAST_VERSION}"
    echo "Download LAST-v${LAST_VERSION}"
    download $LAST_DOWNLOAD_URL "last-${LAST_VERSION}.zip"
    unzip "last-${LAST_VERSION}.zip"
    cd "$build_dir/last-${LAST_VERSION}/src"
    make -j $MAKE_JOBS
    cd ..
    mkdir bin
    cp ./src/lastal ./bin/
    cp ./src/lastal8 ./bin/
    cp ./src/lastdb ./bin/
    cp ./src/lastdb8 ./bin/
    cp ./src/last-merge-batches ./bin/
    cp ./src/last-pair-probs ./bin/
    cp ./src/last-split ./bin/
    cp ./src/last-split8 ./bin/
    cp ./scripts/* ./bin/
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

# --------------- pbmm2 -----------------
echo ""
echo "[$(timestamp)] Installing pbmm2 ..."
pbmm2_dir="$build_dir/pbmm2_conda_env/bin"
if [ -z $(check_installed $pbmm2_dir) ]; then
    cd $build_dir
    clean "$build_dir/pbmm2_conda_env"
    echo "Download pbmm2-v${PBMM2_VERSION}"
    $miniconda3_dir/conda create -y -p $build_dir/pbmm2_conda_env python=3.7
    source $miniconda3_dir/activate $build_dir/pbmm2_conda_env
    $miniconda3_dir/conda install -y -c bioconda/label/cf201901 pbmm2
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
    make modules  
    make -j $MAKE_JOBS
    note_installed $graphmap_dir
fi

# # --------------- GraphMap2 ------------------
# echo ""
# echo "[$(timestamp)] Installing GraphMap2 ..."
# graphmap2_dir="$build_dir/graphmap2-${GRAPHMAP2_VERSION}"
# if [ -z $(check_installed $graphmap2_dir) ]; then
# cd $build_dir
# clean "$build_dir/graphmap2-${GRAPHMAP2_VERSION}"  
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
echo "[$(timestamp)] Installing samtools ..."
samtools_dir="$build_dir/samtools-${SAMTOOLS_VERSION}"
htslib_dir="$samtools_dir/htslib-${SAMTOOLS_VERSION}"
tabix_dir="$samtools_dir/htslib-${SAMTOOLS_VERSION}"
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
    cd htslib-${SAMTOOLS_VERSION}
    ./configure
    make -j $MAKE_JOBS
    cd $build_dir
    rm samtools-${SAMTOOLS_VERSION}.tar.bz2
    note_installed $samtools_dir
fi
PATH="$samtools_dir:$htslib_dir:$tabix_dir:${PATH}"

# --------------- Picard -----------------
echo ""
echo "[$(timestamp)] Installing picard ..."
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
    mkdir GATK3
    cd GATK3
    wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/${SRA_VERSION}/GenomeAnalysisTK.jar
    chmod 755 GenomeAnalysisTK.jar
    ln -s GenomeAnalysisTK.jar gatk3.jar
    note_installed $gatk3_dir
fi

# --------------- GEM-Tools -----------------
echo ""
echo "[$(timestamp)] Installing gemtools ..."                                                                                       
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
echo "[$(timestamp)] Installing freebayes ..."
freebayes_dir="$build_dir/freebayes/bin"
vcflib_dir="$build_dir/freebayes/vcflib/bin"
if [ -z $(check_installed $freebayes_dir) ]; then
    cd $build_dir
    clean "$build_dir/freebayes"
    echo "Download Freebayes-v${FREEBAYES_VERSION}"
    clone $FREEBAYES_DOWNLOAD_URL
    # git clone $FREEBAYES_DOWNLOAD_URL
    cd $build_dir/freebayes
    git checkout -f -q $FREEBAYES_GITHUB_COMMIT_VERSION
    git clone --recursive https://github.com/walaj/SeqLib.git
    cd SeqLib
    git checkout -f -q "5941c68"
    cd ..
    git clone --recursive https://github.com/ekg/bamtools.git
    cd bamtools
    git checkout -f -q "e77a43f"
    cd ..
    git clone --recursive https://github.com/ekg/intervaltree.git
    cd intervaltree
    git checkout -f -q "dbb4c51"
    cd ..
    cd test 
    git clone --recursive https://github.com/illusori/bash-tap.git
    cd bash-tap
    git checkout -f -q "c38fbfa"
    cd ..
    git clone --recursive https://github.com/ingydotnet/test-simple-bash.git
    cd test-simple-bash
    git checkout -f -q "124673f"
    cd ..
    cd ..
    git clone https://github.com/vcflib/vcflib.git
    cd vcflib
    git checkout -f -q "5e3ce04" 
    git clone https://github.com/ekg/fastahack.git
    cd fastahack
    git checkout -f -q "c68cebb"
    cd ..
    git clone https://github.com/ekg/filevercmp.git
    cd filevercmp
    git checkout -f -q "1a9b779"
    cd ..
    git clone https://github.com/ekg/fsom.git
    cd fsom
    git checkout -f -q "a6ef318"
    cd ..
    git clone https://github.com/google/googletest.git
    cd googletest
    git checkout -f -q "d225acc"
    cd ..
    git clone https://github.com/ekg/intervaltree.git
    cd intervaltree
    git checkout -f -q "dbb4c51"
    cd ..
    git clone https://github.com/ekg/multichoose.git
    cd multichoose
    git checkout -f -q "73d35da"
    cd ..
    git clone https://github.com/ekg/smithwaterman.git
    cd smithwaterman
    git checkout -f -q "84c08d7"
    cd ..
    git clone https://github.com/ekg/tabixpp.git
    cd tabixpp
    git checkout -f -q "80012f8"
    git clone https://github.com/samtools/htslib.git
    cd htslib
    git checkout -f -q "0f298ce"
    cd ..
    cd ..
    cd ..
    make 
    cd $build_dir/freebayes/vcflib
    make 
    note_installed $freebayes_dir
fi

# --------------- Clair -----------------
echo ""
echo "[$(timestamp)] Installing Clair ..."
clair_dir="$build_dir/clair_conda_env/bin"
if [ -z $(check_installed $clair_dir) ]; then
    cd $build_dir
    clean "$build_dir/clair_conda_env"
    echo "Download Clair-v${CLAIR_VERSION}"
    $miniconda3_dir/conda create -y -p $build_dir/clair_conda_env python=3.7
    source $miniconda3_dir/activate $build_dir/clair_conda_env
    $miniconda3_dir/conda install -y -c bioconda clair=${CLAIR_VERSION}
    cd $build_dir/clair_conda_env
    # $miniconda3_dir/pypy3 -m ensurepip
    # $miniconda3_dir/pypy3 -m pip install --no-cache-dir intervaltree==3.0.2
    $miniconda3_dir/conda install -y -c conda-forge intervaltree==3.0.2  
    mkdir trained_models && cd trained_models
    $miniconda3_dir/mkdir ont && cd ont
    wget -c --no-check-certificate http://www.bio8.cs.hku.hk/clair_models/ont/122HD34.tar
    tar -xf 122HD34.tar
    cd ..
    mkdir pacbio && cd pacbio
    wget -c --no-check-certificate http://www.bio8.cs.hku.hk/clair_models/pacbio/clr/1234567.tar
    tar -xf 1234567.tar
    cd ..
    mkdir illumina && cd illumina
    wget -c --no-check-certificate http://www.bio8.cs.hku.hk/clair_models/illumina/12345.tar
    tar -xf 12345.tar
    cd ..
    source $miniconda3_dir/deactivate
    note_installed $clair_dir
fi

# --------------- longshot -----------------
echo ""
echo "[$(timestamp)] Installing longshot ..."
longshot_dir="$build_dir/longshot_conda_env/bin"
if [ -z $(check_installed $longshot_dir) ]; then
    cd $build_dir
    clean "$build_dir/longshot_conda_env"
    echo "Download longshot-v${LONGSHOT_VERSION}"
    $miniconda3_dir/conda create -y -p $build_dir/longshot_conda_env python=3.6
    source $miniconda3_dir/activate $build_dir/longshot_conda_env
    $miniconda3_dir/conda install -y -c bioconda longshot=v${LONGSHOT_VERSION}
    source $miniconda3_dir/deactivate
    note_installed $longshot_dir    
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
    $miniconda3_dir/conda create -y -p $build_dir/svim_conda_env python=3.6
    source $miniconda3_dir/activate $build_dir/svim_conda_env
    $miniconda3_dir/conda install -y -c bioconda svim=${SVIM_VERSION}
    source $miniconda3_dir/deactivate
    note_installed $svim_dir
fi

# ------------------ Sniffles ---------------------
echo ""
echo "[$(timestamp)] Installing Sniffles ..."
sniffles_dir="$build_dir/Sniffles-${SNIFFLES_VERSION}/bin/sniffles-core-1.0.11"
if [ -z $(check_installed $sniffles_dir) ]; then
    cd $build_dir
    clean "$build_dir/Sniffles-${SNIFFLES_VERSION}"
    echo "Download Sniffles-v${SNIFFLES_VERSION}"
    download $SNIFFLES_DOWNLOAD_URL "Sniffles-${SNIFFLES_VERSION}.tar.gz"
    tar xvzf Sniffles-${SNIFFLES_VERSION}.tar.gz
    cd $build_dir/Sniffles-${SNIFFLES_VERSION}
    mkdir -p build
    cd build
    cmake ..
    make -j $MAKE_JOBS
    cd $build_dir
    rm Sniffles-${SNIFFLES_VERSION}.tar.gz
    note_installed $sniffles_dir
fi

# ------------------ PBSV ---------------------
echo ""
echo "[$(timestamp)] Installing PBSV ..."
pbsv_dir="$build_dir/pbsv_conda_env/bin"
if [ -z $(check_installed $pbsv_dir) ]; then
    cd $build_dir
    clean "$build_dir/pbsv_conda_env"
    $miniconda3_dir/conda create -y -p $build_dir/pbsv_conda_env python=3.7
    source $miniconda3_dir/activate $build_dir/pbsv_conda_env
    $miniconda3_dir/conda install -y -c bioconda pbsv=${PBSV_VERSION}
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
echo "[$(timestamp)] Installing NanoSV..."
nanosv_dir="$build_dir/nanosv_conda_env/bin"
if [ -z $(check_installed $nanosv_dir) ]; then
    cd $build_dir
    clean "$build_dir/nanosv_conda_env"
    $miniconda3_dir/conda create -y -p $build_dir/nanosv_conda_env python=3.6
    source $miniconda3_dir/activate $build_dir/nanosv_conda_env
    $miniconda3_dir/conda install -y -c bioconda nanosv=${NANOSV_VERSION}
    source $miniconda3_dir/deactivate
    note_installed $nanosv_dir
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

# --------------- ncbi-rmblastn ------------------
echo ""
echo "[$(timestamp)] Installing ncbi-rmblastn ..."
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

# --------------- parallel ------------------
echo ""
echo "[$(timestamp)] Installing parallel ..."
parallel_dir="$build_dir/parallel-${PARALLEL_VERSION}/bin"
if [ -z $(check_installed $parallel_dir) ]; then
    cd $build_dir
    clean "$build_dir/parallel-${PARALLEL_VERSION}"
    echo "Download parallel-${PARALLEL_VERSION}"
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
# Configure executable paths

cd $VARATHON_HOME
echo ""
echo "Configuring executable paths ..."

echo "export VARATHON_HOME=${VARATHON_HOME}" > env.sh
echo "export build_dir=${build_dir}" >> env.sh
echo "export PERL5LIB=${PERL5LIB}" >> env.sh 
echo "export R_LIBS=${R_LIBS}" >> env.sh
echo "export cpanm_dir=${cpanm_dir}" >> env.sh
echo "export miniconda2_dir=${miniconda2_dir}" >> env.sh
echo "export miniconda3_dir=${miniconda3_dir}" >> env.sh
echo "export simuG_dir=${simuG_dir}" >> env.sh
echo "export sra_dir=${sra_dir}" >> env.sh
echo "export art_dir=${art_dir}" >> env.sh
echo "export simlord_dir=${simlord_dir}" >> env.sh
echo "export deepsimulator_dir=${deepsimulator_dir}" >> env.sh
echo "export trimmomatic_dir=${trimmomatic_dir}" >> env.sh
echo "export porechop_dir=${porechop_dir}" >> env.sh
echo "export filtlong_dir=${filtlong_dir}" >> env.sh
echo "export bwa_dir=${bwa_dir}" >> env.sh
echo "export last_dir=${last_dir}" >> env.sh
echo "export ngmlr_dir=${ngmlr_dir}" >> env.sh
echo "export minimap2_dir=${minimap2_dir}" >> env.sh
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
echo "export freebayes_dir=${freebayes_dir}" >> env.sh
echo "export vcflib_dir=${vcflib_dir}" >> env.sh
echo "export vt_dir=${vt_dir}" >> env.sh
echo "export clair_dir=${clair_dir}" >> env.sh
echo "export longshot_dir=${longshot_dir}" >> env.sh
echo "export freec_dir=${freec_dir}" >> env.sh
echo "export manta_dir=${manta_dir}" >> env.sh
echo "export delly_dir=${delly_dir}" >> env.sh
echo "export svim_dir=${svim_dir}" >> env.sh
echo "export sniffles_dir=${sniffles_dir}" >> env.sh
echo "export pbsv_dir=${pbsv_dir}" >> env.sh
echo "export picky_dir=${picky_dir}" >> env.sh
echo "export nanosv_dir=${nanosv_dir}" >> env.sh
echo "export bcftools_dir=${bcftools_dir}" >> env.sh
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
