#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Env;
use Cwd;

##############################################################
#  script: batch_annotate_SNP_INDEL_by_vep.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2022.04.28
#  description: run ensembl-vep annotation for SNVs and INDELs
#  example: perl batch_annotate_SNP_INDEL_by_vep.pl -i Master_Sample_Table.txt -threads 4 -b $batch_id -ref_genome ref_genome.fa  -short_read_mapping_dir ./../01.Short_Read_Mapping -min_mapping_quality 30 -min_variant_calling_quality 30 -caller GATK4  -ploidy 1 -excluded_chr_list yeast.excluded_chr_list.txt
##############################################################

my $VARATHON_HOME = $ENV{VARATHON_HOME};
my $htslib_dir = $ENV{htslib_dir};
my $parallel_dir = $ENV{parallel_dir};
my $vep_dir = $ENV{vep_dir};
my $sample_table = "Master_Sample_Table.txt";
my $batch_id;
my $threads = 1;
my $variant_calling_dir;
my $db_dir;
my $db_cache_version;
my $db_species_name;
my $db_assembly_version;
my $predictor = "vep";
my $mode = "offline";
my $debug = "no";
my $excluded_chr_list;

GetOptions('sample_table|i:s' => \$sample_table,
	   'threads|t:i' => \$threads,
	   'batch_id|b:s' => \$batch_id,
	   'variant_calling_dir|variant_calling_dir:s' => \$variant_calling_dir,
	   'db_dir|db_dir:s' => \$db_dir,
	   'db_cache_version|db_cache_version:s' => \$db_cache_version,
	   'db_species_name|db_species_name:s' => \$db_species_name,
	   'db_assembly_version|db_assembly_version:s' => \$db_assembly_version,
	   'predictor|p:s' => \$predictor,
	   'mode|m:s' => \$mode,
	   'debug|d:s' => \$debug,
	   'excluded_chr_list|excluded_chr_list:s' => \$excluded_chr_list); 


print "Check specified db for variant effect prediction:\n";

if (-e "$db_dir/$db_species_name/${db_cache_version}_${db_assembly_version}") {
    print "Successfully located the specified database: $db_dir/$db_species_name/${db_cache_version}_${db_assembly_version}\n";
} else {
    print "Cannot find the specified database: $db_dir/$db_species_name/${db_cache_version}_${db_assembly_version}\n";
    print "Exit!\n";
    exit;
}

my $sample_table_fh = read_file($sample_table);
my %sample_table = ();
my @sample_table = ();
parse_sample_table($sample_table_fh, \%sample_table, \@sample_table);
my $base_dir = cwd();
my $output_dir = "$batch_id";
system("mkdir $output_dir");


my $sample_count = 0;
foreach my $sample_id (@sample_table) {
    $sample_count++;
    my $local_time = localtime();
    print "[$local_time] processing sample $sample_id \n";
    my $sample_output_dir = "$base_dir/$output_dir/$sample_id";
    system("mkdir -p  $sample_output_dir");
    chdir("$sample_output_dir") or die "cannot change directory to: $!\n";
    system("mkdir tmp");

    print "Check the specified vcf file:\n";

    $local_time = localtime();
    print "[$local_time] processing sample $sample_id for variant effect prediction\n";

    if (-e "$base_dir/$variant_calling_dir/$batch_id/$sample_id") {
        print "Successfully located the specified directory for vcf files: $base_dir/$variant_calling_dir/$batch_id/$sample_id\n";
    } else {
        print "Cannot find the specified directory for vcf files: $base_dir/$variant_calling_dir/$batch_id/$sample_id\n";
        print "Exit!\n";
        exit;
    }

    my @input_vcf = glob "$base_dir/$variant_calling_dir/$batch_id/$sample_id/$sample_id.*.filtered*.vcf";
    foreach my $input_vcf (@input_vcf) {
	my ($vcf_prefix) = ($input_vcf =~ /$base_dir\/$variant_calling_dir\/$batch_id\/$sample_id\/(\S+)\.vcf$/);
	if ($predictor eq "vep") {
	    if ($mode eq "offline") {
		system("/usr/bin/time -v $vep_dir/vep -i $input_vcf --cache --dir_cache $db_dir --species $db_species_name  --assembly $db_assembly_version -o $vcf_prefix.vep.vcf --stats_file $vcf_prefix.vep.summary.html  --force_overwrite --no_check_variants_order --fork $threads --everything --offline");
	    } else {
		system("/usr/bin/time -v $vep_dir/vep -i $input_vcf --cache --dir_cache $db_dir --species $db_species_name  --assembly $db_assembly_version -o $vcf_prefix.vep.vcf --stats_file $vcf_prefix.vep.summary.html  --force_overwrite --no_check_variants_order --fork $threads --everything");
	    }
	} else {
	    die "Unknown predictor: $predictor! Please only use \"vep\" or \"haplosaurus!\n";
	}
    }
    $local_time = localtime();
    print "[$local_time] finishing variant effect prediction for sample $sample_id \n";
	
    # remove intermediate files
    if ($debug eq "no") {

	# system("rm $sample_id.$caller.annotated.vcf");
	# system("rm $sample_id.$caller.annotated.SNP.vcf");
	# system("rm $sample_id.$caller.annotated.INDEL.vcf");
	# system("rm $sample_id.$caller.filtered.SNP.vcf");
	# system("rm $sample_id.$caller.filtered.INDEL.vcf");
	system("rm -r tmp");
    }
    chdir("./../") or die "cannot change directory to: $!\n";
}


print "A total of $sample_count samples were processed!\n";


sub read_file {
    my $file = shift @_;
    my $fh;
    if ($file =~ /\.gz$/) {
        open($fh, "gunzip -c $file |") or die "can't open pipe to $file";
    }
    else {
        open($fh, $file) or die "can't open $file";
    }
    return $fh;
}

sub write_file {
    my $file = shift @_;
    my $fh;
    if ($file =~ /\.gz$/) {
        open($fh, "| gzip -c >$file") or die "can't open $file\n";
    } else {
        open($fh, ">$file") or die "can't open $file\n";
    }
    return $fh;
}  

sub parse_sample_table {
    my ($fh, $sample_table_hashref, $sample_table_arrayref) = @_;
    while (<$fh>) {
        chomp;
        /^#/ and next;
	/^\s*$/ and next;
        my ($sample_id, $PE_read_files, $note) = split /\s+/, $_;
	push @$sample_table_arrayref, $sample_id;
        my ($R1_read_file, $R2_read_file) = split /,/, $PE_read_files;
        $$sample_table_hashref{$sample_id}{'R1_read_file'} = $R1_read_file;
        $$sample_table_hashref{$sample_id}{'R2_read_file'} = $R2_read_file;
	$$sample_table_hashref{$sample_id}{'note'} = $note;
    }
}

