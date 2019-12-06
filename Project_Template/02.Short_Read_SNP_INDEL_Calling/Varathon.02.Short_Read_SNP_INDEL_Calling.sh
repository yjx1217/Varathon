#!/bin/bash
set -e -o pipefail

##########################################
# load environment variables for Varathon
source ./../../env.sh
PATH="$samtools_dir:$PATH"
source $miniconda2_dir/activate $build_dir/conda_clairvoyante_env

###########################################
# set project-specific variables
master_sample_table="Master_Sample_Table.Batch_yeast_sim_illumina.txt" # The master sample table for the processing batch. Default = "Master_Sample_Table.Batch_yeast_sim_illumina.txt".
batch_id="Batch_yeast_sim_illumina" # The batch_id of the processing batch. Default = "Batch_yeast_sim_illumina".
ref_genome="./../00.Reference_Genomes/yeast.tidy.lite.fa" # The input ref_genome used for short-read-based CNV calling. Default = "./../00.Reference_Genomes/yeast.tidy.lite.fa".
excluded_chr_list="./../../data/yeast.excluded_chr_list.txt" # The list for specifying chromosomes/scaffolds/contigs to be exclued for short-read-based CNV calling. Use "" if there is no chromosome/scaffold/contig for exclusion. Default = "./../../data/yeast.excluded_chr_list.txt". 
min_mapping_quality=30 # The minimal mapping quality to use for filtering short-read mapping alignment. Default = 30.
min_variant_calling_quality=30 # The minimal variant calling quality to use for filtering short-read mapping alignment. Default = 30.
ploidy=2; # The ploidy status of samples in the processing batch. Default = 2.
snp_indel_caller="gatk4" # The specific caller used for SNP & INDEL calling: "gatk4" or "freebayes". Default = "gatk4".
threads=4 # The number of threads to use. Defualt = 4.
debug="no" # Whether to keep intermediate files for debuging. Use "yes" if prefer to keep intermediate files, otherwise use "no". Default = "no".
###########################################
# process the pipeline

echo ""
echo "Running short-read-based SNP & INDEL calling ..."
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

short_read_mapping_dir="./../01.Short_Read_Mapping" # The directory for Short-read Mapping. Default = "./../01.Short_Read_Mapping".
perl $VARATHON_HOME/scripts/batch_short_read_SNP_INDEL_calling.pl \
    -i $master_sample_table \
    -t $threads \
    -b $batch_id \
    -ref_genome $ref_genome \
    -short_read_mapping_dir $short_read_mapping_dir \
    -excluded_chr_list $excluded_chr_list \
    -min_mapping_quality $min_mapping_quality \
    -min_variant_calling_quality $min_variant_calling_quality \
    -caller $snp_indel_caller \
    -ploidy $ploidy \
    -debug $debug

echo "Short-read-based SNP & INDEL calling has finished ..."
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
