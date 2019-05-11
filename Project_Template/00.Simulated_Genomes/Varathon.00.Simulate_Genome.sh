#!/bin/bash
set -e -o pipefail

#######################################
# load environment variables for Varathon
source ./../../env.sh

#######################################
# set project-specific variables

ref_genome="./../00.Reference_Genomes/yeast.tidy.lite.fa" # The path to the defined reference genome. Default = "./../00.Reference_Genomes/yeast.tidy.lite.fa" 
excluded_chr_list="./../../data/yeast.excluded_chr_list.txt" # The single-column list defining chromosomes/scaffolds/contigs to be excluded. Default = ./../../data/yeast.excluded_chr_list.txt".
output_prefix="yeast.CNV_DEL" # The filename prefix for the output files. Default = "yeast.CNV_DEL" 
random_seed="20190518" # The random seed used for the genome simulator: simuG.
additional_simuG_options="-cnv_count 10;-cnv_gain_loss_ratio 0;-centromere_gff ./../../data/yeast.centromere.gff3" # For specifying extra options for simuG (https://github.com/yjx1217/simuG.git). Please use semicolons to separate different options.

#######################################
# process the pipeline
echo ""
simuG_options="-refseq $ref_genome;-excluded_chr_list $excluded_chr_list;-prefix $output_prefix;-seed $random_seed;$additional_simuG_options"
echo "Running simuG with the following options:"
OLD_IFS=$IFS
IFS=';'
for op in $simuG_options
do
  echo "$op"
done
IFS=$OLD_IFS

echo ""
echo ""
simuG_options=${simuG_options//;/ }
perl $simuG_dir/simuG.pl $simuG_options


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
