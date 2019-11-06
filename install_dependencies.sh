#!/bin/bash
# last update: 2019.11.05
set -e -o pipefail

VARATHON_HOME=$(pwd)
BUILD="build"

# genome simulation
SIMUG_VERSION="1.0.0" # 
SIMUG_GITHUB_COMMIT_VERSION="472f20c" # commited on 2019.07.12
SIMUG_DOWNLOAD_URL="https://github.com/yjx1217/simuG.git"

# reads retrieving/simulation/processing
SRA_VERSION="2.9.6" # released on 2019.03.18
SRA_DOWNLOAD_URL="https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/${SRA_VERSION}/sratoolkit.${SRA_VERSION}-centos_linux64.tar.gz"

ART_VERSION="mountrainier2016.06.05" # released on 2016.06.05
ART_DOWNLOAD_URL="https://www.niehs.nih.gov/research/resources/assets/docs/artbin${ART_VERSION}linux64.tgz"

SIMLORD_VERSION="1.0.2" #

DEEPSIMULATOR_VERSION="1.5.0" # released on 2019.06.16 
DEEPSIMULATOR_GITHUB_COMMIT_VERSION="3c867c2" # committed on 2019.11.02
DEEPSIMULATOR_DOWNLOAD_URL="https://github.com/lykaust15/DeepSimulator.git"

TRIMMOMATIC_VERSION="0.38" # released on 
TRIMMOMATIC_DOWNLOAD_URL="http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-${TRIMMOMATIC_VERSION}.zip"

PORECHOP_VERSION="0.2.4" # 
PORECHOP_GITHUB_COMMIT_VERSION="109e437" # committed on 2018.10.19
PORECHOP_DOWNLOAD_URL="https://github.com/rrwick/Porechop.git"

FILTLONG_VERSION="0.2.0" #
FILTLONG_GITHUB_COMMIT_VERSION="d1bb46d" # committed on 2018.05.11
FILTLONG_DOWNLOAD_URL="https://github.com/rrwick/Filtlong.git"

# alignment building/processing
BWA_VERSION="0.7.17" # released on 2017.10.23
BWA_DOWNLOAD_URL="https://github.com/lh3/bwa/releases/download/v${BWA_VERSION}/bwa-${BWA_VERSION}.tar.bz2"

LAST_VERSION="979" # released on 
LAST_DOWNLOAD_URL="http://last.cbrc.jp/last-${LAST_VERSION}.zip"

NGMLR_VERSION="0.2.7" # 
NGMLR_DOWNLOAD_URL="https://github.com/philres/ngmlr/releases/download/v${NGMLR_VERSION}/ngmlr-${NGMLR_VERSION}-linux-x86_64.tar.gz"

MINIMAP2_VERSION="2.16" # released on 2019.02.28
MINIMAP2_DOWNLOAD_URL="https://github.com/lh3/minimap2/releases/download/v${MINIMAP2_VERSION}/minimap2-${MINIMAP2_VERSION}_x64-linux.tar.bz2"

PBMM2_VERSION="1.0.0" # released on 2019.03.11
# distributed via conda

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
FREEBAYES_DOWNLOAD_URL="https://github.com/ekg/freebayes.git"

CLAIRVOYANTE_VERSION="1.02" # released on 2019.02.28 

LONGSHOT_VERSION="0.3.4" # released on 2019.05.02

FREEC_VERSION="11.4" #
FREEC_DOWNLOAD_URL="https://github.com/BoevaLab/FREEC/archive/v${FREEC_VERSION}.tar.gz"

MANTA_VERSION="1.5.0" # released on 2018.11.12
MANTA_DOWNLOAD_URL="https://github.com/Illumina/manta/releases/download/v${MANTA_VERSION}/manta-${MANTA_VERSION}.centos6_x86_64.tar.bz2"

DELLY_VERSION="0.8.1" # released on 2019.02.04
DELLY_DOWNLOAD_URL="https://github.com/dellytools/delly/releases/download/v${DELLY_VERSION}/delly_v${DELLY_VERSION}_linux_x86_64bit"

SVIM_VERSION="1.0.0" # 
SVIM_GITHUB_COMMIT_VERSION="b72e631" # committed on 2019.05.02

SNIFFLES_VERSION="1.0.11" #
SNIFFLES_GITHUB_COMMIT_VERSION="e5d5150" # committed on 2019.03.30
SNIFFLES_DOWNLOAD_URL="https://github.com/fritzsedlazeck/Sniffles/archive/${SNIFFLES_VERSION}.tar.gz"

PBSV_VERSION="2.2.0" # released on 2019.02.28

PICKY_VERSION="0.2.a"
PICKY_GITHUB_COMMIT_VERSION="34b85ac" # committed on 2018.07.17
PICKY_DOWNLOAD_URL="https://github.com/TheJacksonLaboratory/Picky.git"

NANOSV_VERSION="1.2.3"
NANOSV_GITHUB_COMMIT_VERSION="c1ae30c" # committed on 2019.04.09
NANOSV_DOWNLOAD_URL="https://github.com/mroosmalen/nanosv.git"

# variant processing
VT_VERSION="" # we use the github commit version below
VT_GITHUB_COMMIT_VERSION="f6d2b5d" # committed on 2018.05.09
VT_DOWNLOAD_URL="https://github.com/atks/vt.git"

HTSLIB_VERSION="1.9" # released on 2018.07.18
HTSLIB_DOWNLOAD_URL="https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2"

BCFTOOLS_VERSION="1.9" # released on 2018.07.18
BCFTOOLS_DOWNLOAD_URL="https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2"

# others
MINICONDA2_VERSION="4.5.11" # released on 2018.09.04
MINICONDA2_DOWNLOAD_URL="https://repo.continuum.io/miniconda/Miniconda2-${MINICONDA2_VERSION}-Linux-x86_64.sh"

BEDTOOLS_VERSION="2.27.1"
BEDTOOLS_DOWNLOAD_URL="https://github.com/arq5x/bedtools2/releases/download/v${BEDTOOLS_VERSION}/bedtools-${BEDTOOLS_VERSION}.tar.gz"

BLAST_VERSION="2.2.31"
BLAST_DOWNLOAD_URL="ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/${BLAST_VERSION}/ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz"

RMBLAST_VERSION="2.2.28"
RMBLAST_DOWNLOAD_URL="ftp://ftp.ncbi.nlm.nih.gov/blast/executables/rmblast/${RMBLAST_VERSION}/ncbi-rmblastn-${RMBLAST_VERSION}-x64-linux.tar.gz"

PARALLEL_VERSION="20180722" # released on 2018.07.22
PARALLEL_DOWNLOAD_URL="http://ftp.gnu.org/gnu/parallel/parallel-${PARALLEL_VERSION}.tar.bz2"

# MUMMER4_VERSION="4.0.0beta2" # released on 2017.10.14
# MUMMER4_DOWNLOAD_URL="https://github.com/gmarcais/mummer/releases/download/v${MUMMER_VERSION}/mummer-${MUMMER_VERSION}.tar.gz"

# GNUPLOT_VERSION="4.6.6"
# GNUPLOT_DOWNLOAD_URL="https://sourceforge.net/projects/gnuplot/files/gnuplot/${GNUPLOT_VERSION}/gnuplot-${GNUPLOT_VERSION}.tar.gz"

# Create the $BUILD directory for dependency installation
if [[ -d $BUILD ]]
then
    echo "The old \"$BUILD\" directory need to deleted for this installation to continue!"
    while true; do
        read -p $'Do you confirm to delete this directory? Please answer yes or no.\n' yn
        case $yn in
            [Yy]* ) echo "The answer \"yes\" is received. The old directory deleted and the new installation start now!"; rm -rf build; break;;
            [Nn]* ) echo "The answer \"no\" is received. Quit the installation!"; exit;;
            * ) echo "Please answer yes or no.";;
        esac
    done
fi

echo ""
echo "Create the new $BUILD directory"
mkdir $BUILD
cd $BUILD
build_dir=$(pwd)

# Downloading all the dependencies
echo ""
echo "Download and install all the dependencies"
echo ""

download () {
  url=$1
  download_location=$2
  echo "Downloading $url to $download_location"
  wget --no-check-certificate $url -O $download_location
}

tidy_version () { 
  echo "$1" | awk -F. '{ printf("%d%03d%03d%03d\n", $1,$2,$3,$4); }';
}


# ---------- set Perl & R environment variables -------------

#PYTHONPATH="$build_dir"
PERL5LIB="$build_dir:$PERL5LIB"
PERL5LIB="$build_dir/cpanm/perlmods/lib/perl5:$PERL5LIB"
R_LIBS="$build_dir/R_libs:$R_LIBS"

mkdir -p  $build_dir/cpanm
cpanm_dir=$build_dir/cpanm
cd $cpanm_dir
wget --no-check-certificate -O - https://cpanmin.us/ > cpanm
chmod +x cpanm
mkdir perlmods
$cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed Test::More@1.302086
$cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed Env@1.04
$cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed Statistics::Descriptive@3.0612
$cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed Statistics::Descriptive::Discrete@0.07
$cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed Math::Random@0.72
$cpanm_dir/cpanm -l $cpanm_dir/perlmods --skip-installed Math::Round@0.07

mkdir $build_dir/R_libs
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

# install dependencies

# ------------- Miniconda --------------------
cd $build_dir
download $MINICONDA2_DOWNLOAD_URL "Miniconda2-${MINICONDA2_VERSION}-Linux-x86_64.sh"
bash Miniconda2-${MINICONDA2_VERSION}-Linux-x86_64.sh -b -p $build_dir/miniconda2
miniconda2_dir="$build_dir/miniconda2/bin"
$miniconda2_dir/conda config --add channels defaults
$miniconda2_dir/conda config --add channels bioconda
$miniconda2_dir/conda config --add channels conda-forge
$miniconda2_dir/pip install cython numpy==1.13.1
rm Miniconda2-${MINICONDA2_VERSION}-Linux-x86_64.sh 

# ------------- simuG -------------------
cd $build_dir
echo "Download simuG-v${SIMUG_VERSION}"
git clone $SIMUG_DOWNLOAD_URL
cd simuG
git checkout -f -q $SIMUG_GITHUB_COMMIT_VERSION
simuG_dir="$build_dir/simuG"

# ------------- SRA Toolkit -------------------
cd $build_dir
echo "Download SRAtoolkit-v${SRA_VERSION}"
download $SRA_DOWNLOAD_URL sratoolkit.${SRA_VERSION}-centos_linux64.tar.gz
tar -zxf sratoolkit.${SRA_VERSION}-centos_linux64.tar.gz
sra_dir="$build_dir/sratoolkit.${SRA_VERSION}-centos_linux64/bin"
rm sratoolkit.${SRA_VERSION}-centos_linux64.tar.gz

# ------------- ART -------------------
cd $build_dir
echo "Download ART-v${ART_VERSION}"
download $ART_DOWNLOAD_URL artbin${ART_VERSION}linux64.tgz
tar -zxf artbin${ART_VERSION}linux64.tgz
art_dir="$build_dir/art_bin_MountRainier"
rm artbin${ART_VERSION}linux64.tgz

# -------------- SimLoRD ----------------
$miniconda2_dir/conda create -y -p $build_dir/conda_simlord_env python=3 pip numpy scipy cython
source $miniconda2_dir/activate $build_dir/conda_simlord_env
pip install pysam
pip install dinopy
pip install "simlord==${SIMLORD_VERSION}"
source $miniconda2_dir/deactivate
simlord_dir="$build_dir/conda_simlord_env/bin"

# -------------- DeepSimulator ----------------
cd $build_dir
echo "Download DeepSimulator-v${DEEPSIMULATOR_GITHUB_COMMIT_VERSION}"
git clone $DEEPSIMULATOR_DOWNLOAD_URL
deepsimulator_dir="$build_dir/DeepSimulator"
cd $deepsimulator_dir
git checkout -f -q $DEEPSIMULATOR_GITHUB_COMMIT_VERSION
$miniconda2_dir/conda remove --name tensorflow_cdpm --all -y
$miniconda2_dir/conda create --name tensorflow_cdpm python=2.7 -y
source $miniconda2_dir/activate tensorflow_cdpm
$miniconda2_dir/conda install -y -c anaconda scikit-learn=0.20.3
pip install numpy==1.13.1
pip install tensorflow==1.2.1
pip install tflearn==0.3.2
pip install tqdm==4.19.4
pip install scipy==0.18.1
pip install h5py==2.7.1
#pip install scikit-learn==0.20.3
pip install biopython==1.74
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

# --------------- Trimmomatic -----------------
cd $build_dir
echo "Download Trimmomatic-v${TRIMMOMATIC_VERSION}"
download $TRIMMOMATIC_DOWNLOAD_URL "Trimmomatic-${TRIMMOMATIC_VERSION}.zip"
unzip Trimmomatic-${TRIMMOMATIC_VERSION}.zip
trimmomatic_dir="$build_dir/Trimmomatic-${TRIMMOMATIC_VERSION}"
cd $trimmomatic_dir
chmod 755 trimmomatic-${TRIMMOMATIC_VERSION}.jar
ln -s trimmomatic-${TRIMMOMATIC_VERSION}.jar trimmomatic.jar 
cd $build_dir
rm Trimmomatic-${TRIMMOMATIC_VERSION}.zip

# --------------- Porechop ------------------
cd $build_dir
echo "Download Porechop-v${PORECHOP_VERSION}"
git clone $PORECHOP_DOWNLOAD_URL
cd Porechop
git checkout -f -q $PORECHOP_GITHUB_COMMIT_VERSION
virtualenv -p $(which python3) py3_virtualenv_porechop
source py3_virtualenv_porechop/bin/activate
py3_virtualenv_porechop/bin/python3 ./setup.py install
deactivate
porechop_dir="$build_dir/Porechop/py3_virtualenv_porechop/bin"

# --------------- Filtlong ------------------
cd $build_dir
echo "Download Filtlong-v${FILTLONG_VERSION}"
git clone $FILTLONG_DOWNLOAD_URL
cd Filtlong
git checkout -f -q $FILTLONG_GITHUB_COMMIT_VERSION
make -j
filtlong_dir="$build_dir/Filtlong/bin"

# ------------- BWA -------------------
cd $build_dir
echo "Download BWA-v${BWA_VERSION}"
download $BWA_DOWNLOAD_URL "bwa-${BWA_VERSION}.tar.bz2"
tar -xjf bwa-${BWA_VERSION}.tar.bz2
bwa_dir="$build_dir/bwa-${BWA_VERSION}"
cd $bwa_dir
make
cd $build_dir
rm bwa-${BWA_VERSION}.tar.bz2

# ------------- LAST ------------------- 
cd $build_dir
echo "Download LAST-v${LAST_VERSION}"
download $LAST_DOWNLOAD_URL "last-${LAST_VERSION}.zip"
unzip "last-${LAST_VERSION}.zip"
cd "$build_dir/last-${LAST_VERSION}/src"
make
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
last_dir="$build_dir/last-${LAST_VERSION}/bin"
cd $build_dir
rm last-${LAST_VERSION}.zip

# ------------- NGMLR ------------------- 
cd $build_dir
echo "Download NGMLR-v${NGMLR_VERSION}"
download $NGMLR_DOWNLOAD_URL "ngmlr-${NGMLR_VERSION}.tar.gz"
tar -xzf "ngmlr-${NGMLR_VERSION}.tar.gz"
ngmlr_dir="$build_dir/ngmlr-${NGMLR_VERSION}"
rm ngmlr-${NGMLR_VERSION}.tar.gz

# --------------- minimap2 ------------------
cd $build_dir
echo "Download minimap2-v${MINIMAP2_VERSION}"
download $MINIMAP2_DOWNLOAD_URL "minimap2-${MINIMAP2_VERSION}.tar.bz2"
tar -xjf minimap2-${MINIMAP2_VERSION}.tar.bz2
minimap2_dir="$build_dir/minimap2-${MINIMAP2_VERSION}_x64-linux"
rm minimap2-${MINIMAP2_VERSION}.tar.bz2

# --------------- pbmm2 -----------------
cd $build_dir
echo "Download pbmm2-v${PBMM2_VERSION}"
$miniconda2_dir/conda create -y -p $build_dir/conda_pbmm2_env
source $miniconda2_dir/activate $build_dir/conda_pbmm2_env
$miniconda2_dir/conda install -y -c bioconda/label/cf201901 pbmm2
source $miniconda2_dir/deactivate
pbmm2_dir="$build_dir/conda_pbmm2_env/bin"

# --------------- samtools -----------------
cd $build_dir
echo "Download samtools-v${SAMTOOLS_VERSION}"
download $SAMTOOLS_DOWNLOAD_URL "samtools-${SAMTOOLS_VERSION}.tar.bz2"
tar -xjf samtools-${SAMTOOLS_VERSION}.tar.bz2
samtools_dir="$build_dir/samtools-${SAMTOOLS_VERSION}"
cd $samtools_dir
C_INCLUDE_PATH=""
./configure --without-curses;
make
cd htslib-${SAMTOOLS_VERSION}
#autoheader
#autoconf
./configure
make
htslib_dir="$samtools_dir/htslib-${SAMTOOLS_VERSION}"
tabix_dir="$samtools_dir/htslib-${SAMTOOLS_VERSION}"
cd $build_dir
rm samtools-${SAMTOOLS_VERSION}.tar.bz2

# --------------- Picard -----------------
cd $build_dir
echo "Download Picard-v${PICARD_VERSION}"
download $PICARD_DOWNLOAD_URL "picard.jar"
mkdir Picard-v${PICARD_VERSION}
picard_dir="$build_dir/Picard-v${PICARD_VERSION}"
mv picard.jar $picard_dir
cd $picard_dir
chmod 755 picard.jar

# --------------- GATK3 ------------------
cd $build_dir
echo "Create GATK3 folder for users' manual installation"
mkdir GATK3
cd GATK3
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/${SRA_VERSION}/GenomeAnalysisTK.jar
chmod 755 GenomeAnalysisTK.jar
ln -s GenomeAnalysisTK.jar gatk3.jar
gatk3_dir="$build_dir/GATK3"

# --------------- GEM-Tools -----------------                                                                                        
cd $build_dir
echo "Download GEMTOOLS-v${GEMTOOLS_VERSION}"
download $GEMTOOLS_DOWNLOAD_URL "gemtools-${GEMTOOLS_VERSION}.tar.gz"
tar zxf "gemtools-${GEMTOOLS_VERSION}.tar.gz"
gemtools_dir="$build_dir/gemtools-${GEMTOOLS_VERSION}-i3/bin"
rm gemtools-${GEMTOOLS_VERSION}.tar.gz

# --------------- GATK4 ------------------                                                                                                
cd $build_dir
echo "Download GATK4-v${GATK4_VERSION}"
download $GATK4_DOWNLOAD_URL "gatk-${GATK4_VERSION}.zip"
unzip gatk-${GATK4_VERSION}.zip
mv gatk-${GATK4_VERSION} GATK4
cd GATK4
ln -s gatk-package-${GATK4_VERSION}-local.jar gatk4.jar
gatk4_dir="$build_dir/GATK4"
cd ..
rm gatk-${GATK4_VERSION}.zip

# --------------- Freebayes -----------------
cd $build_dir
echo "Download Freebayes-v${FREEBAYES_VERSION}"
freebayes_dir="$build_dir/freebayes"
git clone --recursive $FREEBAYES_DOWNLOAD_URL
cd $freebayes_dir
git checkout -f -q $FREEBAYES_GITHUB_COMMIT_VERSION
make
freebayes_dir="$build_dir/freebayes/bin"
vcflib_dir="$build_dir/freebayes/vcflib"
cd $vcflib_dir
make
vcflib_dir="$build_dir/freebayes/vcflib/bin"

# --------------- Clairvoyante -----------------   
cd $build_dir
echo "Download Clairvoyante-v${CLAIRVOYANTE_VERSION}"
$miniconda2_dir/conda create -p $build_dir/conda_clairvoyante_env python=2.7 -y
source $miniconda2_dir/activate $build_dir/conda_clairvoyante_env
$miniconda2_dir/conda install -y -c bioconda clairvoyante=${CLAIRVOYANTE_VERSION}
cd $build_dir/conda_clairvoyante_env
wget https://bootstrap.pypa.io/get-pip.py
pypy get-pip.py
pypy -m pip install --no-cache-dir intervaltree==2.1.0
curl http://www.bio8.cs.hku.hk/trainedModels.tbz | tar -jxf -
source $miniconda2_dir/deactivate
clairvoyante_dir="$build_dir/conda_clairvoyante_env/bin"

# --------------- longshot -----------------   
cd $build_dir
echo "Download longshot-v${LONGSHOT_VERSION}"
$miniconda2_dir/conda create -y -p $build_dir/conda_longshot_env 
source $miniconda2_dir/activate $build_dir/conda_longshot_env
$miniconda2_dir/conda install -y -c bioconda longshot=v${LONGSHOT_VERSION}
source $miniconda2_dir/deactivate
longshot_dir="$build_dir/conda_longshot_env/bin"

# ----------------- FREEC --------------------
cd $build_dir
echo "Download FREEC-v${FREEC_VERSION}"
download $FREEC_DOWNLOAD_URL "FREEC-${FREEC_VERSION}.tar.gz"
tar -xf FREEC-${FREEC_VERSION}.tar.gz
freec_dir="$build_dir/FREEC-${FREEC_VERSION}"
cd $freec_dir/src
make
cp freec $freec_dir
cd $freec_dir/scripts
chmod 755 *
cp * $freec_dir
cd $build_dir
rm FREEC-${FREEC_VERSION}.tar.gz

# ------------------ Manta ---------------------
cd $build_dir
echo "Download Manta-v${MANTA_VERSION}"
download $MANTA_DOWNLOAD_URL "manta-${MANTA_VERSION}.centos6_x86_64.tar.bz2"
tar xjf manta-${MANTA_VERSION}.centos6_x86_64.tar.bz2
manta_dir="$build_dir/manta-${MANTA_VERSION}.centos6_x86_64/bin"
rm manta-${MANTA_VERSION}.centos6_x86_64.tar.bz2

# ------------------ Delly ---------------------
cd $build_dir
echo "Download Delly-v${DELLY_VERSION}"
mkdir delly-${DELLY_VERSION}
cd delly-${DELLY_VERSION}
download $DELLY_DOWNLOAD_URL "delly_v${DELLY_VERSION}_parallel_linux_x86_64bit"
chmod 755 delly_v${DELLY_VERSION}_parallel_linux_x86_64bit
ln -s delly_v${DELLY_VERSION}_parallel_linux_x86_64bit delly
delly_dir="$build_dir/delly-${DELLY_VERSION}"

# ------------------ SVIM ---------------------
cd $build_dir
echo "Download SVIM-v${SVIM_VERSION}"
$miniconda2_dir/conda create -y -p $build_dir/conda_svim_env 
source $miniconda2_dir/activate $build_dir/conda_svim_env
$miniconda2_dir/conda install -y -c bioconda svim=${SVIM_VERSION}
source $miniconda2_dir/deactivate
svim_dir="$build_dir/conda_svim_env/bin"

# ------------------ Sniffles ---------------------
cd $build_dir
echo "Download Sniffles-v${SNIFFLES_VERSION}"
download $SNIFFLES_DOWNLOAD_URL "Sniffles-${SNIFFLES_VERSION}.tar.gz"
tar xzf Sniffles-${SNIFFLES_VERSION}.tar.gz
cd $build_dir/Sniffles-${SNIFFLES_VERSION}
mkdir -p build
cd build
cmake ..
make
sniffles_dir="$build_dir/Sniffles-${SNIFFLES_VERSION}/bin/sniffles-core-1.0.11"
cd $build_dir
rm Sniffles-${SNIFFLES_VERSION}.tar.gz

# ------------------ PBSV ---------------------
cd $build_dir
$miniconda2_dir/conda create -y -p $build_dir/conda_pbsv_env
source $miniconda2_dir/activate $build_dir/conda_pbsv_env
$miniconda2_dir/conda install -y -c bioconda pbsv=${PBSV_VERSION}
source $miniconda2_dir/deactivate
pbsv_dir="$build_dir/conda_pbsv_env/bin"

# ------------------ PICKY ---------------------   
cd $build_dir
echo "Download PICKY-v${PICKY_VERSION}"
git clone $PICKY_DOWNLOAD_URL
cd $build_dir/Picky
git checkout -f -q $PICKY_GITHUB_COMMIT_VERSION
picky_dir="$build_dir/Picky/src"

# ------------------ NANOSV ---------------------
cd $build_dir
$miniconda2_dir/conda create -y -p $build_dir/conda_nanosv_env
source $miniconda2_dir/activate $build_dir/conda_nanosv_env
$miniconda2_dir/conda install -y -c bioconda nanosv=${NANOSV_VERSION}
source $miniconda2_dir/deactivate
nanosv_dir="$build_dir/conda_nanosv_env/bin"

# ------------------ VT ---------------------
cd $build_dir
echo "Download VT-v${VT_VERSION}"
git clone $VT_DOWNLOAD_URL
vt_dir="$build_dir/vt"
cd $vt_dir
git checkout -f -q $VT_GITHUB_COMMIT_VERSION
make
make test

# --------------- bcftools ------------------
cd $build_dir
echo "Download bcftools-v${BCFTOOLS_VERSION}"
download $BCFTOOLS_DOWNLOAD_URL "bcftools-${BCFTOOLS_VERSION}.tar.bz2"
tar -xjf bcftools-${BCFTOOLS_VERSION}.tar.bz2
bcftools_dir="$build_dir/bcftools-${BCFTOOLS_VERSION}"
cd $bcftools_dir
./configure
make
cd $build_dir
rm bcftools-${BCFTOOLS_VERSION}.tar.bz2

# --------------- bedtools ------------------
cd $build_dir
echo "Download bedtools-v${BEDTOOLS_VERSION}"
download $BEDTOOLS_DOWNLOAD_URL "bedtools-${BEDTOOLS_VERSION}.tar.gz"
tar -zxf bedtools-${BEDTOOLS_VERSION}.tar.gz
cd "$build_dir/bedtools2"
make
bedtools_dir="$build_dir/bedtools2/bin"
cd $build_dir
rm bedtools-${BEDTOOLS_VERSION}.tar.gz

# --------------- ncbi-blast+ ------------------
cd $build_dir
echo "Download ncbi-blast-v${BLAST_VERSION}"
download $BLAST_DOWNLOAD_URL "ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz"
tar -zxf ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz
blast_dir="$build_dir/ncbi-blast-${BLAST_VERSION}+/bin"
windowmasker_dir="$build_dir/ncbi-blast-${BLAST_VERSION}+/bin"
rm ncbi-blast-${BLAST_VERSION}+-x64-linux.tar.gz

# --------------- ncbi-rmblast ------------------
cd $build_dir
echo "Download ncbi-rmblastn-v${BLAST_VERSION}"
download $RMBLAST_DOWNLOAD_URL "ncbi-rmblastn-${RMBLAST_VERSION}-x64-linux.tar.gz"
tar -zxf ncbi-rmblastn-${RMBLAST_VERSION}-x64-linux.tar.gz
rmblast_dir="$build_dir/ncbi-rmblastn-${RMBLAST_VERSION}/bin"
# copy rmblastn binary file to ncbi-blast+ directory for easy RepeatMasker configuration
cp $rmblast_dir/rmblastn $blast_dir
rm ncbi-rmblastn-${RMBLAST_VERSION}-x64-linux.tar.gz

# --------------- parallel ------------------
parallel_dir="$build_dir/parallel-${PARALLEL_VERSION}/bin"
cd $build_dir
echo "Download parallel"
download $PARALLEL_DOWNLOAD_URL "parallel_v${PARALLEL_VERSION}.tar.bz2"
tar -xvjf parallel_v${PARALLEL_VERSION}.tar.bz2
cd parallel-${PARALLEL_VERSION}
./configure --prefix="$build_dir/parallel-${PARALLEL_VERSION}"
make
make install
parallel_dir="$build_dir/parallel-${PARALLEL_VERSION}/bin"
cd ..
rm parallel_v${PARALLEL_VERSION}.tar.bz2

#############
# # --------------- mummer4 ------------------
# cd $build_dir
# echo "Download mummer-v${MUMMER_VERSION}"
# download $MUMMER_DOWNLOAD_URL "mummer-${MUMMER_VERSION}.tar.gz"
# tar -zxf mummer-${MUMMER_VERSION}.tar.gz
# mummer4_dir="$build_dir/mummer-${MUMMER_VERSION}"
# cd $mummer4_dir
# ./configure
# make
# PATH=$mummer4_dir:$PATH
# cd $build_dir
# rm mummer-${MUMMER_VERSION}.tar.gz

# # --------------- gnuplot ------------------
# cd $build_dir
# echo "Download gnuplot-v${GNUPLOT_VERSION}"
# download $GNUPLOT_DOWNLOAD_URL "gnuplot-${GNUPLOT_VERSION}.tar.gz"
# tar -zxf gnuplot-${GNUPLOT_VERSION}.tar.gz
# cd "$build_dir/gnuplot-${GNUPLOT_VERSION}"
# ./configure --prefix="$build_dir/gnuplot-${GNUPLOT_VERSION}" --disable-wxwidgets
# make
# make install
# gnuplot_dir="$build_dir/gnuplot-${GNUPLOT_VERSION}/bin"
# cd $build_dir
# rm gnuplot-${GNUPLOT_VERSION}.tar.gz
# PATH="$gnuplot_dir:${PATH}"

# Configure executable paths

cd $VARATHON_HOME

echo "Configuring executable paths ..."

echo "export VARATHON_HOME=${VARATHON_HOME}" > env.sh
echo "export build_dir=${build_dir}" >> env.sh
echo "export PERL5LIB=${PERL5LIB}" >> env.sh 
echo "export R_LIBS=${R_LIBS}" >> env.sh
echo "export cpanm_dir=${cpanm_dir}" >> env.sh
echo "export miniconda2_dir=${miniconda2_dir}" >> env.sh
echo "export simuG_dir=${simuG_dir}" >> env.sh
echo "export sra_dir=${sra_dir}" >> env.sh
echo "export art_dir=${art_dir}" >> env.sh
echo "export simlord_dir=${simlord_dir}" >> env.sh
echo "export nanosim_dir=${nanosim_dir}" >> env.sh
echo "export deepsimulator_dir=${deepsimulator_dir}" >> env.sh
echo "export trimmomatic_dir=${trimmomatic_dir}" >> env.sh
echo "export porechop_dir=${porechop_dir}" >> env.sh
echo "export filtlong_dir=${filtlong_dir}" >> env.sh
echo "export bwa_dir=${bwa_dir}" >> env.sh
echo "export last_dir=${last_dir}" >> env.sh
echo "export ngmlr_dir=${ngmlr_dir}" >> env.sh
echo "export minimap2_dir=${minimap2_dir}" >> env.sh
echo "export pbmm2_dir=${pbmm2_dir}" >> env.sh
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
echo "export clairvoyante_dir=${clairvoyante_dir}" >> env.sh
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
    echo ""
    echo "Varathon message: This bash script has been successfully processed! :)"
    echo ""
    echo ""
    exit 0
fi
############################
