#!/bin/bash
set -e -o pipefail

##########################################
# load environment variables for Varathon
source ./../../env.sh

###########################################
# set project-specific variables
master_sample_table="Master_Sample_Table.Batch_yeast_sim_pacbio.txt" # The master sample table for the processing batch. Default = "Master_Sample_Table.Batch_yeast_sim_pacbio.txt".
batch_id="Batch_yeast_sim_pacbio" # The batch_id of the processing batch. Default = "Batch_yeast_sim_pacbio".
ref_genome="./../00.Reference_Genomes/yeast.tidy.lite.fa" # The input ref_genome used for long-read-based SV calling. Default = "./../00.Reference_Genomes/yeast.tidy.lite.fa".
excluded_chr_list="./../../data/yeast.excluded_chr_list.txt" # The list for specifying chromosomes/scaffolds/contigs to be exclued for long-read-based CNV calling. Use "" if there is no chromosome/scaffold/contig for exclusion. Default = "./../../data/yeast.excluded_chr_list.txt". 
min_mapping_quality=30 # The minimal mapping quality to use for filtering long-read mapping alignment. Default = 30.
min_read_support=5 # The minimal number of reads supporting the corresponding SV. Default = 5. 
long_read_technology="pacbio" # The sequencing technology of simulated long reads. E.g. "pacbio" or "nanopore". Default = "pacbio".
sv_caller="sniffles" # The specific caller used for SV calling. Supported callers include: "sniffles", "svim", "picky", "pbsv", and "nanosv". Default = "sniffles".
threads=4 # The number of threads to use. Defualt = 4.
debug="no" # Whether to keep intermediate files for debuging. Use "yes" if prefer to keep intermediate files, otherwise use "no". Default = "no".
###########################################
# process the pipeline

echo ""
echo "Running long-read-based SV calling ..."
echo ""
test_file_existence () {
    filename=$1
    if [[ ! -f $filename ]]
    then
        echo "The file $filename does not exists! Process terminated!"
        exit
    else
        echo "Test passed!"
    fi
}

echo "Testing the existence of ref_genome: $ref_genome .."
test_file_existence $ref_genome

echo ""
long_read_mapping_dir="./../11.Long_Read_Mapping" # The directory for Long-read Mapping. Default = "./../11.Long_Read_Mapping".
# ploidy=1 # The ploidy status of samples in the processing batch. Default = 2.

perl $VARATHON_HOME/scripts/batch_long_read_SV_calling.pl \
    -i $master_sample_table \
    -t $threads \
    -b $batch_id \
    -ref_genome $ref_genome \
    -long_read_mapping_dir $long_read_mapping_dir \
    -excluded_chr_list $excluded_chr_list \
    -min_mapping_quality $min_mapping_quality \
    -min_read_support $min_read_support \
    -long_read_technology $long_read_technology \
    -caller $sv_caller \
    -debug $debug

echo "Long-read-based SV calling has finished ..."
echo ""
# clean up intermediate files
echo ""
echo "Clean up intermediate files ..."
echo ""
if [[ $debug = "no" ]]
then
    :
fi

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
