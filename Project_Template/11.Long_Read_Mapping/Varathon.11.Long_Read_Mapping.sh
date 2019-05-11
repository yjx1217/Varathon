#!/bin/bash
set -e -o pipefail

##########################################
# load environment variables for Varathon
source ./../../env.sh

###########################################
# set project-specific variables

master_sample_table="Master_Sample_Table.Batch_yeast_sim_pacbio.txt" # The master sample table for the processing batch. Default = "Master_Sample_Table.Batch_yeast_sim_pacbio.txt".
batch_id="Batch_yeast_sim_pacbio" # The batch_id used for the processing batch. Default = "Batch_yeast_sim_pacbio".
ref_genome="./../00.Reference_Genomes/yeast.tidy.lite.fa" # The input reference genome to be used for short-read mapping. Default = "./../00.Reference_Genomes/yeast.tidy.lite.fa".
excluded_chr_list="./../../data/yeast.excluded_chr_list.txt" # The list for specifying chromosomes/scaffolds/contigs to be exclued for long-read mapping. Use "" if there is no chromosome/scaffold/contig for exclusion. Default = "./../../data/yeast.excluded_chr_list.txt". 
long_read_dir="./../00.Long_Reads" # The directory for the input long reads. Default = "./../00.Long_Reads".
long_read_technology="pacbio" # The long-read technology used: "pacbio" or "nanopore". Default = "pacbio".
long_read_mapper="minimap2" # The long-read mapper used for mapping. Supported mappers include: "minimap2", "ngmlr", "last", and "pbmm2". Default = "minimap2".
min_mapping_quality=30 # The minimal mapping quality to use for filtering long-read mapping alignment. Default = 30.
threads=4 # The number of threads to use. Defualt = 4.
debug="no" # Whether to keep intermediate files for debuging. Use "yes" if prefer to keep intermediate files, otherwise use "no". Default = "no".
###########################################
# process the pipeline

echo ""
echo "Running long-read mapping for $long_read_technology reads with $long_read_mapper ..."
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


perl $VARATHON_HOME/scripts/batch_long_read_mapping.pl \
    -i $master_sample_table \
    -t $threads \
    -ref_genome $ref_genome \
    -long_read_dir $long_read_dir \
    -long_read_mapper $long_read_mapper \
    -long_read_technology $long_read_technology \
    -batch_id $batch_id \
    -excluded_chr_list $excluded_chr_list \
    -min_mapping_quality $min_mapping_quality \
    -debug $debug

echo "Long-read mapping has finished ..."
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
