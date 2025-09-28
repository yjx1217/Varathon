#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Env;
use Cwd;

##############################################################
#  script: batch_long_read_mapping.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2020.05.28
#  description: run long-read mapping for a batch of samples
#  example: perl batch_long_read_mapping.pl -i Master_Sample_Table.txt -b $batch_id -threads 4  -ref_genome ref_genome.fa -long_read_dir ./../00.Long_Reads -long_read_technology pacbio -long_read_mapper minimap2 -min_mapping_quality 30 -excluded_chr_list yeast.excluded_chr_list.txt
##############################################################

my $VARATHON_HOME = $ENV{VARATHON_HOME};
my $java_dir = $ENV{java_dir};
my $minimap2_dir = $ENV{minimap2_dir};
my $ngmlr_dir = $ENV{ngmlr_dir};
my $last_dir = $ENV{last_dir};
my $picky_dir = $ENV{picky_dir};
my $pbmm2_dir = $ENV{pbmm2_dir};
my $graphmap_dir = $ENV{graphmap_dir};
my $graphmap2_dir = $ENV{graphmap2_dir};
my $samtools_dir = $ENV{samtools_dir};
my $picard_dir = $ENV{picard_dir};
my $sample_table = "Master_Sample_Table.txt";
my $threads = 1;
my $batch_id = "Batch_TEST";
my $long_read_technology = "pacbio";
my $long_read_mapper = "minimap2";
my $min_mapping_quality = 20;
my $ref_genome;
my $long_read_dir;
my $excluded_chr_list;
my $debug = "no";

GetOptions('sample_table|i:s' => \$sample_table,
	   'threads|t:i' => \$threads,
	   'batch_id|b:s' => \$batch_id,
	   'ref_genome|ref_genome:s' => \$ref_genome,
	   'long_read_dir|long_read_dir:s' => \$long_read_dir,
	   'long_read_technology|long_read_technology:s' => \$long_read_technology,
	   'long_read_mapper|m:s' => \$long_read_mapper,
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

print "Check the specified reference genome file:\n";
my $ref_genome_file = "$base_dir/$ref_genome";
if (-e $ref_genome_file) {
    print "Successfully located the specified reference genome file: $ref_genome_file.\n";
} else {
    print "Cannot find the specified reference genome file: $ref_genome_file!\n";
    print "Exit!\n";
    exit;
}

if (-d "$base_dir/$output_dir/ref_genome_preprocessing") {
    print "found pre-calculated ref_geome_preprocessing directory!\n";
} else {
    print "generating ref_geome_preprocessing directory!\n";
    system("mkdir -p $base_dir/$output_dir/ref_genome_preprocessing");
    chdir("$base_dir/$output_dir/ref_genome_preprocessing");
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
}



foreach my $sample_id (@sample_table) {
    my $local_time = localtime();
    print "[$local_time] processing sample $sample_id with long-read mapping\n";
    $sample_count++;

    my $long_read_file = "$base_dir/$long_read_dir/$sample_table{$sample_id}{'long_read_file'}";
    print "Check the specified long read file:\n";
    if (-e $long_read_file) {
	print "Successfully located the specified long read file: $long_read_file.\n";
    } else {
	print "Cannot find the specified long read file: $long_read_file!\n";
	print "Exit!\n";
	exit;
    }
    print "processing sample $sample_id\n";
    print "long_read_file = $long_read_file\n";
    my $sample_output_dir;
    print "mapping $long_read_file to reference genome: $ref_genome\n";
    $mapping_count++;
    $sample_output_dir = "$base_dir/$output_dir/$sample_id";
    system("mkdir -p  $sample_output_dir");
    chdir("$sample_output_dir") or die "cannot change directory to: $!\n";
    system("mkdir tmp");
    print("mapping the reads by $long_read_mapper\n");
    system("cp $base_dir/$output_dir/ref_genome_preprocessing/ref.genome.fa .");
    system("cp $base_dir/$output_dir/ref_genome_preprocessing/ref.genome.fa.fai .");
    system("cp $base_dir/$output_dir/ref_genome_preprocessing/ref.genome.dict .");
    if ($long_read_mapper eq "minimap2") {
	if ($long_read_technology eq "pacbio") {
	    system("/usr/bin/time -v $minimap2_dir/minimap2 -t $threads -ax map-pb ref.genome.fa $base_dir/$long_read_dir/$sample_table{$sample_id}{'long_read_file'} | $samtools_dir/samtools view -bS -q $min_mapping_quality - >$sample_id.bam");
	} elsif ($long_read_technology eq "nanopore") {
	    system("/usr/bin/time -v $minimap2_dir/minimap2 -t $threads -ax map-ont ref.genome.fa $base_dir/$long_read_dir/$sample_table{$sample_id}{'long_read_file'} | $samtools_dir/samtools view -bS -q $min_mapping_quality - >$sample_id.bam");
	} else {
	    system("/usr/bin/time -v $minimap2_dir/minimap2 -t $threads -ax map-hifi ref.genome.fa $base_dir/$long_read_dir/$sample_table{$sample_id}{'long_read_file'} | $samtools_dir/samtools view -bS -q $min_mapping_quality - >$sample_id.bam");
	}
	if ($debug eq "no") {
	    # system("rm ref.genome.fa.amb");
	    # system("rm ref.genome.fa.sa");
	}
    } elsif ($long_read_mapper eq "ngmlr") {
	if ($long_read_technology eq "pacbio") {
	    system("/usr/bin/time -v $ngmlr_dir/ngmlr -t $threads -x pacbio -r ref.genome.fa -q $base_dir/$long_read_dir/$sample_table{$sample_id}{'long_read_file'} | $samtools_dir/samtools view -bS -q $min_mapping_quality - >$sample_id.bam");
	} else {
	    system("/usr/bin/time -v $ngmlr_dir/ngmlr -t $threads -x ont -r ref.genome.fa -q $base_dir/$long_read_dir/$sample_table{$sample_id}{'long_read_file'} | $samtools_dir/samtools view -bS -q $min_mapping_quality - >$sample_id.bam");
	}
	if ($debug eq "no") {
	    system("rm ref.genome.fa-*.ngm");
	}
    } elsif ($long_read_mapper eq "last") {
	system("/usr/bin/time -v $last_dir/lastdb -cR01 -v -P $threads ref.genome.lastdb ref.genome.fa");
	system("/usr/bin/time -v $last_dir/lastal -C 2 -K 2 -r 1 -q 3 -a 2 -b 1 -Q 1 -i 100M -P $threads ref.genome.lastdb $base_dir/$long_read_dir/$sample_table{$sample_id}{'long_read_file'} |gzip -c > $sample_id.maf.gz");
	# system("/usr/bin/time -v $last_dir/lastal -r 1 -q 1 -a 0 -b 2 -Q 1 -i 100M -P $threads ref.genome.lastdb $base_dir/$long_read_dir/$sample_table{$sample_id}{'long_read_file'} |gzip -c > $sample_id.maf.gz");
	system("/usr/bin/time -v $java_dir/java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar CreateSequenceDictionary -REFERENCE ref.genome.fa -OUTPUT ref.genome.dict");
	system("/usr/bin/time -v gzip -dc $sample_id.maf.gz | $last_dir/last-split - |gzip -c > $sample_id.lastsplit.maf.gz");
	system("/usr/bin/time -v gzip -dc $sample_id.lastsplit.maf.gz | $last_dir/maf-convert -f ref.genome.dict sam -r \"ID:$sample_id PL:$long_read_technology SM:$sample_id\" -  | $samtools_dir/samtools view -bS -q $min_mapping_quality - >$sample_id.bam");
        if ($debug eq "no") {
	    system("rm ref.genome.lastdb.tis");
	    system("rm ref.genome.lastdb.ssp");
	    system("rm ref.genome.lastdb.sds");
	    system("rm ref.genome.lastdb.des");
	    system("rm ref.genome.lastdb.suf");
	    system("rm ref.genome.lastdb.prj");
	    system("rm ref.genome.lastdb.bck");
        }
    } elsif ($long_read_mapper eq "pbmm2") {
	system("gzip -c -d $base_dir/$long_read_dir/$sample_table{$sample_id}{'long_read_file'} > $sample_id.raw.fastq"); 
	system("/usr/bin/time -v $pbmm2_dir/pbmm2 align --sort --preset SUBREAD --sample $sample_id -j $threads ref.genome.fa $sample_id.raw.fastq $sample_id.bam");
	if ($debug eq "no") {
	    system("rm $sample_id.raw.fastq");
        }
    } elsif ($long_read_mapper eq "graphmap") {
	system("gzip -c -d $base_dir/$long_read_dir/$sample_table{$sample_id}{'long_read_file'} > $sample_id.raw.fastq");
	system("/usr/bin/time -v $graphmap_dir/graphmap align -t $threads -r ref.genome.fa -d $sample_id.raw.fastq -o $sample_id.sam");
        if ($debug eq "no") {
            system("rm $sample_id.raw.fastq");
            system("rm ref.genome.fa.gmidx");
        }
    } else {
	die "Error! Unrecognized long_read_mapper: $long_read_mapper! Exit! \n";
    }
    if ($long_read_mapper =~ /(minimap2|ngmlr|pbmm2|last|graphmap)/) {
	if ($long_read_mapper =~ /(minimap2|ngmlr|last|graphmap)/) {
	    ## sort bam file by picard-tools: SortSam
	    system("$java_dir/java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar SortSam -INPUT $sample_id.bam -OUTPUT $sample_id.sort.bam -SORT_ORDER coordinate -VALIDATION_STRINGENCY LENIENT -MAX_RECORDS_IN_RAM 50000");
	} else {
	    system("$samtools_dir/samtools view -bS -q $min_mapping_quality $sample_id.bam >$sample_id.sort.bam");
	}
	if ($debug eq "no") {
	    if (-e "$sample_id.bam") {
		system("rm $sample_id.bam");
	    }
	    if (-e "$sample_id.bam.bai") {
		system("rm $sample_id.bam.bai");
	    }
	}
	# Picard tools remove duplicates (majorly for amplicon sequencing)
	system("$java_dir/java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar MarkDuplicates -INPUT $sample_id.sort.bam -REMOVE_DUPLICATES true -METRICS_FILE $sample_id.dedup.matrics -VALIDATION_STRINGENCY LENIENT -OUTPUT $sample_id.dedup.bam"); 
	if ($debug eq "no") {
	    system("rm $sample_id.sort.bam");
	    system("rm $sample_id.dedup.matrics");
	}
	# index final bam file
	system("$samtools_dir/samtools index $sample_id.dedup.bam");
	# generate samtools mpileup 
	# system("$samtools_dir/samtools mpileup -Q 0 -C 50 -q $min_mapping_quality -f ref.genome.fa $sample_id.dedup.bam |gzip -c >${sample_id}.mpileup.gz");	
	# compute basic alignment statistics by samtools
	system("$samtools_dir/samtools flagstat $sample_id.dedup.bam >$sample_id.samstat");
	# calculate per-base depth
	system("$samtools_dir/samtools depth -aa $sample_id.dedup.bam |gzip -c >$sample_id.depth.txt.gz");
	# calculate read mapping coverage statistics
	system("perl $VARATHON_HOME/scripts/summarize_mapping_coverage.pl -r ref.genome.fa -s $sample_id.samstat -d $sample_id.depth.txt.gz -c 5 -t $sample_id -o $sample_id.coverage_summary.txt");
	if ($mapping_count == 1) {
	    system("cp $sample_id.coverage_summary.txt ./../all_samples.coverage_summary.txt");
	} else {
	    system("tail -1 $sample_id.coverage_summary.txt  >> ./../all_samples.coverage_summary.txt");
	}
	system("rm -r tmp");
	$local_time = localtime();
	print "[$local_time] finishing long-read mapping for sample $sample_id \n";
  
	# remove intermediate files
	if ($debug eq "no") {
	    system("rm ref.genome.fa");
	    system("rm ref.genome.fa.fai");
	    system("rm ref.genome.dict");
	    # system("rm $sample_id.depth.txt.gz");
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
        my ($sample_id, $long_read_file, $note) = split /\s+/, $_;
	push @$sample_table_arrayref, $sample_id;
	$$sample_table_hashref{$sample_id}{'long_read_file'} = $long_read_file;
	$$sample_table_hashref{$sample_id}{'note'} = $note;
    }
}

