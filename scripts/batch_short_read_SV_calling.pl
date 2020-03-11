#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Env;
use Cwd;

##############################################################
#  script: batch_short_read_SV_calling.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2020.04.27
#  description: run short-read-based SV calling for a batch of samples
#  example: perl batch_short_read_SV_calling.pl -i Master_Sample_Table.txt -threads 4 -b $batch_id -ref_genome ref_genome.fa  -short_read_mapping_dir ./../01.Short_Read_Mapping -min_mapping_quality 20 -caller manta  -ploidy 1 -excluded_chr_list yeast.excluded_chr_list.txt
##############################################################

my $VARATHON_HOME = $ENV{VARATHON_HOME};
my $java_dir = $ENV{java_dir};
my $samtools_dir = $ENV{samtools_dir};
my $picard_dir = $ENV{picard_dir};
my $manta_dir = $ENV{manta_dir};
my $delly_dir = $ENV{delly_dir};
my $bcftools_dir = $ENV{bcftools_dir};
my $tabix_dir = $ENV{tabix_dir};
my $sample_table = "Master_Sample_Table.txt";
my $batch_id;
my $threads = 1;
my $ref_genome;
my $short_read_mapping_dir;
my $min_mapping_quality = 30;
my $min_variant_calling_quality = 30;
my $ploidy = 2;
my $caller = "manta"; # "manta" or "delly"
my $debug = "no";
my $excluded_chr_list;

GetOptions('sample_table|i:s' => \$sample_table,
	   'threads|t:i' => \$threads,
	   'batch_id|b:s' => \$batch_id,
	   'ref_genome|ref_genome:s' => \$ref_genome,
	   'short_read_mapping_dir|short_read_mapping_dir:s' => \$short_read_mapping_dir,
	   'min_mapping_quality|mq:i' => \$min_mapping_quality,
	   'caller|c:s' => \$caller,
	   'ploidy|p:i' => \$ploidy,
	   'debug|d:s' => \$debug,
	   'excluded_chr_list|excluded_chr_list:s' => \$excluded_chr_list); 

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
    print "[$local_time] processing sample $sample_id with read-alignment filtering\n";
    my $ref_genome_file = "$base_dir/$ref_genome";
    my $bam_file = "$base_dir/$short_read_mapping_dir/$batch_id/$sample_id/$sample_id.realn.bam";
    print "Check the specified reference genome file:\n";
    if (-e $ref_genome_file) {
        print "Successfully located the specified reference genome file: $ref_genome_file.\n";
    } else {
        print "Cannot find the specified reference genome file: $ref_genome_file!\n";
        print "Exit!\n";
        exit;
    }
    print "Check the specified short-read mapping bam file:\n";
    if (-e $bam_file) {
        print "Successfully located the specified short-read mapping bam file: $bam_file.\n";
    } else {
        print "Cannot find the specified short-read mapping bam file: $bam_file!\n";
        print "Exit!\n";
        exit;
    }
    my $sample_output_dir = "$base_dir/$output_dir/$sample_id";
    system("mkdir -p  $sample_output_dir");
    chdir("$sample_output_dir") or die "cannot change directory to: $!\n";
    system("mkdir tmp");
    if (defined $excluded_chr_list) {
	if (-e "$base_dir/$excluded_chr_list") {
	    system("$VARATHON_HOME/scripts/select_fasta_by_list.pl -i $base_dir/$ref_genome -l $base_dir/$excluded_chr_list -m reverse -o ref.genome.fa");
	} else {
	    die "cannot find $excluded_chr_list at $base_dir/$excluded_chr_list\n";
	}
    } else {
        system("cp $base_dir/$ref_genome ref.genome.fa");
    }
    ## index reference sequence
    system("$samtools_dir/samtools faidx ref.genome.fa");
    system("$java_dir/java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar CreateSequenceDictionary -REFERENCE ref.genome.fa -OUTPUT ref.genome.dict");
    ## filter bam file by mapping quality
    system("$samtools_dir/samtools view -b -q $min_mapping_quality $base_dir/$short_read_mapping_dir/$batch_id/$sample_id/$sample_id.realn.bam  >$sample_id.filtered.bam");
    # index the filtered.bam file
    system("$samtools_dir/samtools index $sample_id.filtered.bam");

    $local_time = localtime();
    print "[$local_time] processing sample $sample_id with SV calling\n";

    # SV calling 
    if ($caller eq "manta") {
	system("/usr/bin/time -v $manta_dir/configManta.py --bam $sample_id.filtered.bam --referenceFasta ref.genome.fa --runDir ${sample_id}_manta_out");
	system("/usr/bin/time -v ./${sample_id}_manta_out/runWorkflow.py -m local -j $threads");
	system("cp ./${sample_id}_manta_out/results/variants/diploidSV.vcf.gz $sample_id.manta.SV.vcf.gz");
	system("gunzip $sample_id.manta.SV.vcf.gz");
	#system("cp ./${sample_id}_manta_out/results/variants/diploidSV.vcf.gz.tbi $sample_id.manta.SV.vcf.gz.tbi");
    } elsif ($caller eq "delly") {
	system("/usr/bin/time -v $delly_dir/delly call -g ref.genome.fa -o $sample_id.delly.SV.bcf $sample_id.filtered.bam");
	system("/usr/bin/time -v $bcftools_dir/bcftools view $sample_id.delly.SV.bcf > $sample_id.delly.SV.vcf");
	#system("$tabix_dir/bgzip -c $sample_id.delly.SV.vcf > $sample_id.delly.SV.vcf.gz");
	#system("$tabix_dir/tabix -p vcf $sample_id.delly.SV.vcf.gz");
    } else {
	die "Unknown caller: $caller! Please use either \"manta\" or \"delly\" as your short-read-based SV caller!\n";
    }

    $local_time = localtime();
    print "[$local_time] finishing short-read-based SV calling for sample $sample_id \n";
  
    # remove intermediate files
    if ($debug eq "no") {
	system("rm ref.genome.fa");
	system("rm ref.genome.fa.fai");
	system("rm ref.genome.dict");
	system("rm $sample_id.filtered.bam");
	system("rm $sample_id.filtered.bam.bai");
	system("rm -r tmp");
	if ($caller eq "delly") {
	    system("rm $sample_id.delly.SV.bcf");
	    system("rm $sample_id.delly.SV.bcf.csi");
	}
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

sub parse_sample_cnv_file {
    my ($sample_cnv_fh, $sample, $refseq, $all_samples_cnv_fh) = @_;
    while (<$sample_cnv_fh>) {
	chomp;
	/^#/ and next;
	/^\s*$/ and next;
	/chr\tstart/ and next;
	my ($chr, $start, $end, $copy_number, $status, $MWU_test_p_value, $KS_test_p_value) = split /\t/, $_;
	if ($start <= $end) {
	    print $all_samples_cnv_fh "$sample\t$refseq\t$_\n";
	}
    }
}

sub get_read_length {
    my $fh = shift @_;
    my $read_length = 100;
    while (<$fh>) {
        chomp;
        if ($. % 4 == 2) {
            $read_length = length $_;
            last;
        }
    }
    return $read_length;
}

