#!/bin/bash
set -e -o pipefail

#######################################
# load environment variables for Varathon
source ./../../env.sh

#######################################
# set project-specific variables

ref_genome="./../00.Simulated_Genomes/yeast.CNV_DEL.simseq.genome.fa" # Reference or Simulated genome to be used for simulating reads. Default = ./../00.Simulated_Genomes/yeast.CNV_DEL.simseq.genome.fa".
long_read_technology="pacbio" # The sequencing technology of simulated long reads. E.g. "pacbio" or "nanopore". Default = "pacbio".
read_coverage="30" # The total coverage of simulated reads. Default = "30".
output_prefix="yeast.CNV_DEL" # Filename prefix for the output files. Default = "yeast.CNV_DEL".
threads="4" # The number of threads to use. Default = "4".
random_seed="20190518" # Random seed for read simulation. Default = 20190518".
debug="no" # Wether to keep intermediate files for debugging. Choose "yes" if you want to keep the intermediate files. Otherwise use "no". Default = "no".
#######################################

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

echo "Simulating $long_read_technology reads for the genome $ref_genome .." 
if [[ $long_read_technology == "pacbio" ]]
then
    $simlord_dir/simlord \
	--read-reference $ref_genome \
	--coverage $read_coverage \
	--no-sam \
	$output_prefix.$long_read_technology
else
    export PATH="$VARATHON_HOME/build/miniconda2/bin:$PATH"
    $deepsimulator_dir/deep_simulator.sh -H $deepsimulator_dir -c $threads -i $ref_genome -B 2 -K $read_coverage -S $random_seed -o ${output_prefix}_deepsimulator_out
    cp ./${output_prefix}_deepsimulator_out/pass.fastq $output_prefix.$long_read_technology.fastq
    if [[ $debug == "no" ]]
    then
	rm -r ${output_prefix}_deepsimulator_out
    fi
fi

echo ""
echo "Reads simulation is done!"
echo "Compress the simulated reads"
gzip $output_prefix.$long_read_technology.fastq
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
