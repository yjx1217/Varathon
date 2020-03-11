#!/bin/bash
set -e -o pipefail

#######################################
# load environment variables for Varathon
source ./../../env.sh

#######################################
# set project-specific variables
sra_run_list="./../../data/yeast.long_reads.sra_run_list.txt" # A simple tab separated file with two compulsory columns, in which the first column contains the SRA run id and the sencond column contains the corresponding sample name. Lines started with "#" will be ignored.
#######################################
# process the pipeline

# check if defined sra_run_list exists
test_file_existence () {
    filename=$1
    if [[ ! -f $filename ]]
    then
        echo "The file $filename does not exists! Process terminated!"
        exit
    fi
}

test_file_existence $sra_run_list

# mkdir tmp
i=0
while read -r line
do
    [[ $line == \#* ]] && continue
    [[ $line == "" ]] && continue
    IFS=$'\t' read -r srr_id sample_name <<<"$line"
    i=$((i+1))
    if [[ -z $sample_name ]]
    then
	sample_name="sample_${i}"
    fi
    echo "retrieve reads by the SRR_id: $srr_id for the sample $sample_name ..."
    # $sra_dir/fasterq-dump --threads $threads --temp ./tmp $srr_id 
    $sra_dir/fastq-dump --defline-seq '@$sn[_$rn]/$ri' --defline-qual '+$sn[_$rn]/$ri' \
     	--gzip -skip-technical --dumpbase --read-filter pass --clip $srr_id
    ln -s ${srr_id}_pass.fastq.gz $sample_name.fastq.gz    
done < $sra_run_list
# rm -r tmp

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
