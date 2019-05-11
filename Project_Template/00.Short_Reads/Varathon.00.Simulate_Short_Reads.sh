#!/bin/bash
set -e -o pipefail

#######################################
# load environment variables for Varathon
source ./../../env.sh

#######################################
# set project-specific variables
ref_genome="./../00.Simulated_Genomes/yeast.CNV_DEL.simseq.genome.fa" # Reference or Simulated genome to be used for simulating reads. Default = ./../00.Simulated_Genomes/yeast.CNV_DEL.simseq.genome.fa".
Illumina_platform_and_read_length="HiSeq2500L150" # The Illumina sequencing platform and the associated read length. Available options include: "HiSeq2kL100" for HiSeq2000 platform with 100-bp read length; "HiSeq2500L150" for HiSeq2500 platform with 150-bp read length; "HiSeq2500L150" for HiSeq2500 platform with 125-bp read length. "HiSeqXtruSeqL150" for HiSeqX TrueSeq platform with 150-bp read length; "MiSeqv3L250" for MiSeq platform with 250-bp read length; and "NextSeq500v2L75" for NextSeq platform with 75-bp read length. Default = "HiSeq2500L150".
read_coverage="30" # The coverage of simulated Illumina reads. Default = "30".
output_prefix="yeast.CNV_DEL" # Filename prefix for the output files. Default = "yeast.CNV_DEL".
random_seed="20190518" # Random seed used for read simulation. Default = 20190518".

#######################################

# threads="1" # The number of threads to use. Default = "1".

# process the pipeline
echo ""
test_file_existence () {
    filename=$1
    if [[ ! -f $filename ]]
    then
	echo "The file $filename does not exists! Process terminated!"
	exit
    fi
}

test_file_existence $ref_genome

echo "Simulating Illumina reads for the genome $ref_genome .." 

read1_profile="$art_dir/Illumina_profiles/${Illumina_platform_and_read_length}R1.txt"
read2_profile="$art_dir/Illumina_profiles/${Illumina_platform_and_read_length}R2.txt"

test_file_existence $read1_profile
test_file_existence $read2_profile

echo "Simulating Illumina reads with the following read profile:"
echo "$read1_profile"
echo "$read2_profile"

read_length=${Illumina_platform_and_read_length##*L}
echo "read_length=$read_length" 
echo ""

mean_fragment_size=500

$art_dir/art_illumina \
    --qprof1 $read1_profile \
    --qprof2 $read2_profile \
    -f $read_coverage \
    -l $read_length \
    -i $ref_genome \
    -p \
    -na \
    -rs $random_seed \
    -m $mean_fragment_size \
    -s 10 \
    -o $output_prefix.illumina.R

echo ""
echo "Reads simulation is done!"

mv ${output_prefix}.illumina.R1.fq ${output_prefix}.illumina.R1.fastq
mv ${output_prefix}.illumina.R2.fq ${output_prefix}.illumina.R2.fastq

echo "Compress the simulated reads"

gzip ${output_prefix}.illumina.R1.fastq
gzip ${output_prefix}.illumina.R2.fastq
# $htslib_dir/bgzip --threads $threads < $output_prefix.illumina.R1.fastq > $output_prefix.illumina.R1.fastq.gz
# $htslib_dir/bgzip --threads $threads < $output_prefix.illumina.R2.fastq > $output_prefix.illumina.R2.fastq.gz

echo ""
echo "Done!"
echo ""

############################
# checking bash exit status
if [[ $? -eq 0 ]]
then
    echo ""
    echo "Varathon message: This bash script has been successfully processed! :)"
    echo ""
    echo ""
    exit 0
fi
############################
