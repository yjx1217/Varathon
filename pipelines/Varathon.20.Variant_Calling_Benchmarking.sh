#!/bin/bash
set -e -o pipefail

##########################################
# load environment variables for Varathon
source ./../../env.sh

###########################################
# set project-specific variables
ref_vcf="./../00.Simulated_Genomes/simuG_yeast/yeast.SNP.refseq2simseq.SNP.vcf.gz" # The path to the reference VCF file to be used as the ground truth.
query_vcf="./../02.Short_Read_SNP_INDEL_Calling/Batch_yeast_sim_illumina/yeast.SNP/yeast.SNP.gatk4.filtered.SNP.vcf.gz" # The path to the query VCF file to be used for comparing with the reference VCF file. 
vcf_for_sv="no" # Whether the input VCF files are for denoting structural variants (SVs). Default = "no".

ref_genome="./../00.Reference_Genomes/yeast.tidy.lite.fa" # The reference genome used for SNP/INDEL calling. Needed when vcf_for_sv="no". Default = "./../00.Reference_Genomes/yeast.tidy.lite.fa"
snp_indel_qual_filter=30 # The SNP/INDEL variant quality filter to be used when comparing reference and query VCF files. Needed when vcf_for_sv="no". Default = "30".
sv_qual_filter="LowQual;q5" # The SV variant quality filter to be used when comparing reference and query VCF files. Needed when vcf_for_sv="yes". Default = "LowQual;q5".
max_offset=10 # The maximal breakend coordinate offset allowance to be used when comparing reference and query vcf files. Needed when vcf_for_sv="yes". Default = 10. (i.e. 10 bp).

output_prefix="yeast.SNP" # The filename prefix for the output report. Default = "yeast.SNP".
debug="no" # Whether to keep intermediate files for debuging. Use "yes" if prefer to keep intermediate files, otherwise use "no". Default = "no".
###########################################
# process the pipeline

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

echo "Testing the existence of ref_vcf: $ref_vcf .."
test_file_existence $ref_vcf
echo "Testing the existence of query_vcf: $query_vcf .."
test_file_existence $query_vcf

if [[ $vcf_for_sv == "no" ]]
then
    echo "Testing the existence of ref_genome: $ref_genome .."
    test_file_existence $ref_genome
fi

echo ""
if [[ $ref_vcf =~ \.gz$ ]]
then
    gunzip -c $ref_vcf > ref_vcf.raw.vcf
else
    cp $ref_vcf ref_vcf.raw.vcf
fi

if [[ $query_vcf =~ \.gz$ ]]
then
    gunzip -c $query_vcf > query_vcf.raw.vcf
else
    cp $query_vcf query_vcf.raw.vcf
fi

if [[ $vcf_for_sv == "no" ]]
then
    echo "Normalizing input vcf files ..."
    $vt_dir/vt decompose_blocksub ref_vcf.raw.vcf -a -o ref_vcf.decompose.vcf
    $vt_dir/vt normalize ref_vcf.decompose.vcf -r $ref_genome -f "~VARIANT_CONTAINS_N"| $vt_dir/vt uniq - -o  ref_vcf.normalize.vcf
    $vt_dir/vt decompose_blocksub query_vcf.raw.vcf -a -o query_vcf.decompose.vcf
    $vt_dir/vt normalize query_vcf.decompose.vcf -r $ref_genome -f "~VARIANT_CONTAINS_N"| $vt_dir/vt uniq - -o  query_vcf.normalize.vcf
    echo "######################"
    echo "ref_vcf=$ref_vcf"
    echo "query_vcf=$query_vcf"
    echo "vcf_for_sv=$vcf_for_sv"
    echo "snp_indel_qual_filter=$snp_indel_qual_filter"
    echo ""
    perl $VARATHON_HOME/scripts/vcf_benchmarker.pl \
	-ref_vcf ref_vcf.normalize.vcf \
	-query_vcf query_vcf.normalize.vcf \
	-vcf_for_sv $vcf_for_sv \
	-snp_indel_qual_filter $snp_indel_qual_filter \
	-p $output_prefix
else
    echo "######################"
    echo "ref_vcf=$ref_vcf"
    echo "query_vcf=$query_vcf"
    echo "vcf_for_sv=$vcf_for_sv"
    echo "max_offset=$max_offset"
    echo "sv_qual_filter=$sv_qual_filter"
    echo ""
    perl $VARATHON_HOME/scripts/vcf_benchmarker.pl \
	-ref_vcf ref_vcf.raw.vcf \
	-query_vcf query_vcf.raw.vcf \
	-vcf_for_sv $vcf_for_sv \
	-max_offset $max_offset \
	-sv_qual_filter $sv_qual_filter \
	-p $output_prefix
fi

if [[ $debug == "no" ]]
then
    rm ref_vcf.raw.vcf
    rm query_vcf.raw.vcf
    if [[ $vcf_for_sv == "no" ]]
    then
	rm ref_vcf.decompose.vcf
	rm ref_vcf.normalize.vcf
	rm query_vcf.decompose.vcf
	rm query_vcf.normalize.vcf
    fi
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

