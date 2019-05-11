#!/bin/bash
set -e -o pipefail

##########################################
# load environment variables for Varathon
source ./../../env.sh

###########################################
# set project-specific variables
master_sample_table="Master_Sample_Table.Batch_yeast_sim_illumina.txt" # The master sample table for the processing batch. Default = "Master_Sample_Table.Batch_yeast_sim_illumina.txt"
batch_id="Batch_yeast_sim_illumina" # The batch_id used for the processing batch. Default = "Batch_yeast_sim_illumina".
ref_genome="./../00.Reference_Genomes/yeast.tidy.lite.fa" # The input reference genome to be used for short-read mapping. Default = "./../00.Reference_Genomes/yeast.tidy.lite.fa".
excluded_chr_list="./../../data/yeast.excluded_chr_list.txt" # The list for specifying chromosomes/scaffolds/contigs to be exclued for short-read mapping. Use "" if there is no chromosome/scaffold/contig for exclusion. Default = "./../../data/yeast.excluded_chr_list.txt". 
short_read_dir="./../00.Short_Reads" # The directory for the input short reads. Default = "./../00.Short_Reads".
min_mapping_quality=30 # The minimal mapping quality to use for filtering short-read mapping alignment. Default = 30.
threads=4 # The number of threads to use. Defualt = 4.
debug="no" # Whether to keep intermediate files for debuging. Use "yes" if prefer to keep intermediate files, otherwise use "no". Default = "no".
###########################################
# process the pipeline

echo ""
echo "Running short-read mapping ..."
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

perl $VARATHON_HOME/scripts/batch_short_read_mapping.pl \
    -i $master_sample_table \
    -t $threads \
    -ref_genome $ref_genome \
    -short_read_dir $short_read_dir \
    -batch_id $batch_id \
    -excluded_chr_list $excluded_chr_list \
    -min_mapping_quality $min_mapping_quality \
    -debug $debug

echo "Short-read mapping has finished ..."
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
