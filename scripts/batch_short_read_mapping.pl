#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Env;
use Cwd;

##############################################################
#  script: batch_short_read_mapping.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2020.04.27
#  description: run short-read mapping for a batch of samples
#  example: perl batch_short_read_mapping.pl -i Master_Sample_Table.txt -b $batch_id -threads 4  -ref_genome ref_genome.fa -short_read_dir ./../00.Offspring_Reads -min_mapping_quality 30 -excluded_chr_list yeast.excluded_chr_list.txt
##############################################################

my $VARATHON_HOME = $ENV{VARATHON_HOME};
my $java_dir = $ENV{java_dir};
my $trimmomatic_dir = $ENV{trimmomatic_dir};
my $bwa_dir = $ENV{bwa_dir};
my $samtools_dir = $ENV{samtools_dir};
my $picard_dir = $ENV{picard_dir};
my $gatk3_dir = $ENV{gatk3_dir};
my $sample_table = "Master_Sample_Table.txt";
my $threads = 1;
my $batch_id = "Batch_TEST";
my $min_mapping_quality = 20;
my $ref_genome;
my $short_read_dir;
my $excluded_chr_list;
my $debug = "no";

GetOptions('sample_table|i:s' => \$sample_table,
	   'threads|t:i' => \$threads,
	   'batch_id|b:s' => \$batch_id,
	   'ref_genome|ref_genome:s' => \$ref_genome,
	   'short_read_dir|short_read_dir:s' => \$short_read_dir,
	   'min_mapping_quality|q:i' => \$min_mapping_quality,
	   'debug|d:s' => \$debug,
	   'excluded_chr_list|excluded_chr_list:s' => \$excluded_chr_list); 

my $output_dir = "$batch_id";
my $sample_table_fh = read_file($sample_table);
my %sample_table = ();
my @sample_table = ();
parse_sample_table($sample_table_fh, \%sample_table, \@sample_table);
my $base_dir = cwd();
system("mkdir $output_dir");

my $sample_count = 0;
my $mapping_count = 0;
my $adapter = "$trimmomatic_dir/adapters/TruSeq3-PE-2.fa";

foreach my $sample_id (@sample_table) {
    my $local_time = localtime();
    print "[$local_time] processing sample $sample_id with short-read mapping\n";
    $sample_count++;
    my $ref_genome_file = "$base_dir/$ref_genome";
    my $R1_read_file = "$base_dir/$short_read_dir/$sample_table{$sample_id}{'R1_read_file'}";
    my $R2_read_file = "$base_dir/$short_read_dir/$sample_table{$sample_id}{'R2_read_file'}";
    print "Check the specified reference genome file:\n";
    if (-e $ref_genome_file) {
        print "Successfully located the specified reference genome file: $ref_genome_file.\n";
    } else {
        print "Cannot find the specified reference genome file: $ref_genome_file!\n";
        print "Exit!\n";
        exit;
    }
    print "Check the specified short read file:\n";
    if (-e $R1_read_file) {
        print "Successfully located the specified short read file 1: $R1_read_file.\n";
    } else {
        print "Cannot find the specified short read file 1: $R1_read_file!\n";
        print "Exit!\n";
        exit;
    }
    if (-e $R2_read_file) {
        print "Successfully located the specified short read file 2: $R2_read_file.\n";
    } else {
        print "Cannot find the specified short read file 2: $R2_read_file!\n";
        print "Exit!\n";
        exit;
    }

    print "processing sample $sample_id\n";
    print "PE_read_file = $R1_read_file,$R2_read_file\n";
    my $R1_read_fh = read_file($R1_read_file);
    my $raw_read_length = get_read_length($R1_read_fh);
    print "raw read length = $raw_read_length\n";
    my $sample_output_dir;
    print "mapping $R1_read_file,$R2_read_file to reference genome: $ref_genome\n";
    $mapping_count++;
    $sample_output_dir = "$base_dir/$output_dir/$sample_id";
    system("mkdir -p  $sample_output_dir");
    chdir("$sample_output_dir") or die "cannot change directory to: $!\n";
    system("mkdir tmp");
    print "trim the reads by trimmomatic\n";
    system("ln -s $adapter adapter.fa");
    system("ln -s $R1_read_file $sample_id.R1.raw.fq.gz");
    system("ln -s $R2_read_file $sample_id.R2.raw.fq.gz");
    system("/usr/bin/time -v $java_dir/java -Djava.io.tmpdir=./tmp -XX:ParallelGCThreads=$threads -jar $trimmomatic_dir/trimmomatic.jar PE -threads $threads -phred33 $sample_id.R1.raw.fq.gz  $sample_id.R2.raw.fq.gz $sample_id.R1.trimmed.PE.fq.gz $sample_id.R1.trimmed.SE.fq.gz  $sample_id.R2.trimmed.PE.fq.gz $sample_id.R2.trimmed.SE.fq.gz  ILLUMINACLIP:adapter.fa:2:30:10  SLIDINGWINDOW:5:20 MINLEN:36");
    system("rm $sample_id.R1.raw.fq.gz");
    system("rm $sample_id.R2.raw.fq.gz");
    system("rm $sample_id.R1.trimmed.SE.fq.gz");
    system("rm $sample_id.R2.trimmed.SE.fq.gz");
    system("rm adapter.fa");
    print("mapping the reads by bwa\n");
    if (defined $excluded_chr_list) {
	if (-e "$base_dir/$excluded_chr_list") { 
	    system("$VARATHON_HOME/scripts/select_fasta_by_list.pl -i $base_dir/$ref_genome -l $base_dir/$excluded_chr_list -m reverse -o ref.genome.fa");
	} else {
	    die "$excluded_chr_list not found at $base_dir/$excluded_chr_list\n";
	}
    } else {
	system("cp $base_dir/$ref_genome ref.genome.fa");
    }
    system("$bwa_dir/bwa index ref.genome.fa");
    system("/usr/bin/time -v $bwa_dir/bwa mem -t $threads -M ref.genome.fa $sample_id.R1.trimmed.PE.fq.gz $sample_id.R2.trimmed.PE.fq.gz >$sample_id.sam");
    if ($debug eq "no") {
	system("rm ref.genome.fa.bwt");
	system("rm ref.genome.fa.pac");
	system("rm ref.genome.fa.ann");
	system("rm ref.genome.fa.amb");
	system("rm ref.genome.fa.sa");
	system("rm $sample_id.R1.trimmed.PE.fq.gz");
	system("rm $sample_id.R2.trimmed.PE.fq.gz");
    }
    ## index reference sequence
    system("$samtools_dir/samtools faidx ref.genome.fa");
    system("$java_dir/java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar CreateSequenceDictionary -REFERENCE ref.genome.fa -OUTPUT ref.genome.dict");
    ## filter bam file by mapping quality
    system("$samtools_dir/samtools view -bS -q $min_mapping_quality $sample_id.sam >$sample_id.bam");
    if ($debug eq "no") {
	system("rm $sample_id.sam");
    }
    ## sort bam file by picard-tools: SortSam
    system("$java_dir/java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar SortSam -INPUT $sample_id.bam -OUTPUT $sample_id.sort.bam -SORT_ORDER coordinate");
    if ($debug eq "no") {
	system("rm $sample_id.bam");
    }
    ## fixmate
    system("$java_dir/java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar FixMateInformation -INPUT $sample_id.sort.bam -OUTPUT $sample_id.fixmate.bam");
    if ($debug eq "no") {
	system("rm $sample_id.sort.bam");
    }
    ## add or replace read groups and sort
    system("$java_dir/java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar AddOrReplaceReadGroups -INPUT $sample_id.fixmate.bam -OUTPUT $sample_id.rdgrp.bam -SORT_ORDER coordinate -RGID $sample_id -RGLB $sample_id -RGPL 'Illumina' -RGPU $sample_id -RGSM $sample_id -RGCN 'RGCN'");
    if ($debug eq "no") {
	system("rm $sample_id.fixmate.bam");
    }
    # Picard tools remove duplicates
    system("$java_dir/java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar MarkDuplicates -INPUT $sample_id.rdgrp.bam -REMOVE_DUPLICATES true -METRICS_FILE $sample_id.dedup.matrics -OUTPUT $sample_id.dedup.bam"); 
    if ($debug eq "no") {
	system("rm $sample_id.rdgrp.bam");
    }
    # index the dedup.bam file
    system("$samtools_dir/samtools index $sample_id.dedup.bam");
    # GATK local realign
    # find realigner targets
    system("$java_dir/java -Djava.io.tmpdir=./tmp -XX:ParallelGCThreads=$threads -jar $gatk3_dir/gatk3.jar -nt $threads -R ref.genome.fa -T RealignerTargetCreator -I $sample_id.dedup.bam  -o $sample_id.realn.intervals");
    # run realigner
    system("$java_dir/java -Djava.io.tmpdir=./tmp -XX:ParallelGCThreads=$threads -jar $gatk3_dir/gatk3.jar -R ref.genome.fa -T IndelRealigner -I $sample_id.dedup.bam -targetIntervals $sample_id.realn.intervals -o $sample_id.realn.bam");
    if ($debug eq "no") {
    	system("rm $sample_id.dedup.bam");
    	system("rm $sample_id.dedup.bam.bai");
    	system("rm $sample_id.dedup.matrics");
    }
    # index final the bam file
    system("$samtools_dir/samtools index $sample_id.realn.bam");
    if ($debug eq "no") {
    	system("rm $sample_id.realn.intervals");
    }
    # generate samtools mpileup 
   system("$samtools_dir/samtools mpileup -Q 0 -C 50 -q $min_mapping_quality -f ref.genome.fa  $sample_id.realn.bam |gzip -c >${sample_id}.mpileup.gz");	
    # compute basic alignment statistics by samtools
    system("$samtools_dir/samtools flagstat $sample_id.realn.bam >$sample_id.samstat");
    # calculate per-base depth
    system("$samtools_dir/samtools depth -aa $sample_id.realn.bam |gzip -c >$sample_id.depth.txt.gz");
    # compute insert size statistics
    # system("$java_dir/java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar CollectInsertSizeMetrics -INPUT $sample_id.realn.bam -OUTPUT $sample_id.insert_size_metrics.txt -H $sample_id.insert_size_histogram.pdf -M 0.5");
    # calculate read mapping coverage statistics
    system("perl $VARATHON_HOME/scripts/summarize_mapping_coverage.pl -r ref.genome.fa -s $sample_id.samstat -d $sample_id.depth.txt.gz -c 5 -t $sample_id -o $sample_id.coverage_summary.txt");
    if ($mapping_count == 1) {
	system("cp $sample_id.coverage_summary.txt ./../all_samples.coverage_summary.txt");
    } else {
	system("tail -1 $sample_id.coverage_summary.txt  >> ./../all_samples.coverage_summary.txt");
    }
    system("rm -r tmp");
    $local_time = localtime();
    print "[$local_time] finishing short-read mapping for sample $sample_id \n";
  
    # remove intermediate files
    if ($debug eq "no") {
	system("rm ref.genome.fa");
	system("rm ref.genome.fa.fai");
	system("rm ref.genome.dict");
	# system("rm $sample_id.depth.txt.gz");
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

