# Varathon

<p align="center">
  <img src="https://github.com/yjx1217/Varathon/blob/master/Varathon.logo.png" alt="Varathon logo" width="545" height="200"/>
</p>

**Varathon: a scalable variant calling and benchmarking framework supporting both short and long reads**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Description
Varathon is a scalable variant calling and benchmarking framework that supports both short and long reads. In addition to small variants such as SNPs and INDELs, Varathon can also identify large variants such as structural variants (SVs) (e.g. inversions, translocations, segmental deletions and duplications). Copy-number variants (CNVs), the results of segmental deletions and duplications, can also be nicely profiled using a traditional sliding-window-based method with short-reads. 

Under the hood, a series of task-specific modules are provided to carry out the full workflow of read-mapping-based variant calling:

* **00. Reference_Genomes**
  * preparing reference genome(s)
* **00.Simulated_Genomes**
  * generating simulated genome(s) (by simuG)
* **00.Short_Reads**
  * downloading (by SRA tools) or simulating (by ART) short reads
* **00.Long_Reads**
  * downloading (by SRA tools) or simulating (by SimLoRD or DeepSimulator) long reads
* **01.Short_Read_Mapping**
  * mapping short-read mapping to the reference genome (by bwa)
* **02.Short_Read_SNP_INDEL_Calling**
  * calling SNPs and INDELs based on the short-read mapping alignment (by GATK4, freebayes, or clair)
* **03.Short_Read_SV_Calling**
  * calling SVs based on the short-read mapping alignment (by Manta or Delly)
* **04.Short_Read_CNV_Calling**
  * calling CNVs based on the short-read mapping alignment (by FREEC+DNAcopy)
* **11.Long_Read_Mapping**
  * mapping long reads to the reference genome (by minimap2, ngmlr, last, or pbmm2)
* **12.Long_Read_SNP_INDEL_Calling**
  * calling SNPs and INDELs based on the long-read mapping alignment (by longshot or clair)
* **13.Long_Read_SV_Calling**
  * calling SVs based on the long-read mapping alignment (by Sniffles, svim, Picky, NanoSV, or pbsv)
* **20.Variant_Calling_Benchmarking**
  * comparing different variant calling VCF files to calculate benchmarking statistics such as precision, recall, and F1 score.

## License
Varathon is distributed under the MIT license.

## Installation
```sh
git clone https://github.com/yjx1217/Varathon.git
cd Varathon
bash ./install_dependencies.sh
```

## Requirements
### Hardware, operating system and network requirements
Varathon is designed for a desktop or computing server running an x86-64-bit Linux operating system. Multithreaded processors are preferred to speed up the process since many modules support multithreaded processing. 

### Software requirements
* bash (https://www.gnu.org/software/bash/)
* bzip2 and libbz2-dev (http://www.bzip.org/)
* curl (https://curl.haxx.se/)
* gcc and g++ (https://gcc.gnu.org/)
* git (https://git-scm.com/)
* GNU make (https://www.gnu.org/software/make/)
* gzip (https://www.gnu.org/software/gzip/)
* libopenssl-devel
* libcurl-devel
* java runtime environment (JRE) v1.8.0 (https://www.java.com)
* perl v5.12 or newer (https://www.perl.org/)
* python v2.7.9 or newer (https://www.python.org/)
* python v3.4 or newer (https://www.python.org/)
* tar (https://www.gnu.org/software/tar/)
* unzip (http://infozip.sourceforge.net/UnZip.html)
* wget (https://www.gnu.org/software/wget/)
* zlib and zlib-devel (https://zlib.net/)
* xz and xz-devel (https://tukaani.org/xz/)
