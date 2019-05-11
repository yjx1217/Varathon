#!/bin/bash
set -e -o pipefail

#######################################
# load environment variables for Varathon
source ./../../env.sh

#######################################
# set project-specific variables

ref_genome_prefix="yeast" # The file name prefix of the reference genome. Default = "yeast". 
ref_genome_download_URL="ftp://ftp.ensembl.org/pub/release-95/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa.gz" # The URL for downloading the reference genome. 
excluded_chr_list="./../../data/yeast.excluded_chr_list.txt" # The single-column list defining chromosomes/scaffolds/contigs to be excluded. Use "" if there is no chromosome/scaffold/contig for exclusion. Default = ./../../data/yeast.excluded_chr_list.txt".

#######################################
# process the pipeline

download_and_extract() {
    url=$1
    echo "Downloading $url"
    if [[ $url =~ \.gz$ ]]; 
    #if [[ $url =~ \.fa.gz$ || $url =~ \.fasta.gz$ ]]; 
    then
	download_location="$ref_genome_prefix.raw.fa.gz"
        extract_command="gunzip"
	wget --no-check-certificate $url -O $download_location
	gunzip $download_location
    else
	download_location="$ref_genome_prefix.raw.fa"
	wget --no-check-certificate $url -O $download_location
    fi
}

echo ""
echo "Retrieve the sample reference genome assembly ..."
download_and_extract $ref_genome_download_URL 
echo ""
echo "Tidy the sample reference genome assembly ..."
$VARATHON_HOME/scripts/tidy_fasta.pl -i $ref_genome_prefix.raw.fa -o $ref_genome_prefix.tidy.fa
echo ""
echo "Removing the excluded chromosomes ..."
if [[ ! -z "$excluded_chr_list" ]]
then
    $VARATHON_HOME/scripts/select_fasta_by_list.pl -i $ref_genome_prefix.tidy.fa -l $excluded_chr_list -m reverse -o $ref_genome_prefix.tidy.lite.fa
else
    cp $ref_genome_prefix.tidy.fa $ref_genome_prefix.tidy.lite.fa
fi

echo ""
echo "Removing intermediate files ..."
rm $ref_genome_prefix.raw.fa
rm $ref_genome_prefix.tidy.fa

echo ""
echo "Done!"

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
