#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Statistics::Descriptive;
use Env;
use Cwd;

##############################################################
#  script: batch_short_read_CNV_calling.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2023.05.03
#  description: run short-read-based CNV calling for a batch of samples
#  example: perl batch_short_read_CNV_calling.pl -i Master_Sample_Table.txt -threads 4 -b $batch_id -ref_genome ref_genome.fa  -short_read_mapping_dir ./../01.Short_Read_Mapping -min_mapping_quality 30 -window_size 250  -ploidy 1 -excluded_chr_list yeast.excluded_chr_list.txt
##############################################################

my $VARATHON_HOME = $ENV{VARATHON_HOME};
my $java_dir = $ENV{java_dir};
my $samtools_dir = $ENV{samtools_dir};
my $picard_dir = $ENV{picard_dir};
my $bedtools_dir = $ENV{bedtools_dir};
my $gemtools_dir = $ENV{gemtools_dir};
my $freec_dir = $ENV{freec_dir};
my $sample_table = "Master_Sample_Table.txt";
my $batch_id;
my $threads = 1;
my $min_mapping_quality = 30;
my $ref_genome;
my $short_read_mapping_dir;
my $min_mappability = 0.85;
my $window_size = 250;
my $step_size = $window_size;
my $ploidy = 2;
my $excluded_chr_list;
my $debug = "no";

GetOptions('sample_table|i:s' => \$sample_table,
	   'threads|t:i' => \$threads,
	   'batch_id|b:s' => \$batch_id,
	   'min_mapping_quality|q:i' => \$min_mapping_quality,
	   'ref_genome|ref_genome:s' => \$ref_genome,
	   'short_read_mapping_dir|short_read_mapping_dir:s' => \$short_read_mapping_dir,
	   'min_mappability|min_map:f' => \$min_mappability,
	   'window_size|w:i' => \$window_size,
	   'step_size|s:i' => \$step_size,
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

my $all_samples_cnv = "$base_dir/$output_dir/all_samples.CNV_summary.txt";
my $all_samples_cnv_fh = write_file($all_samples_cnv);
print $all_samples_cnv_fh "sample\tref_genome\tchr\tstart\tend\tcopy_number\tstatus\tMWU_test_p_value\tKS_test_p_value\n";

my $sample_count = 0;
my $raw_read_length = 100; # we fix read length = 100 bp for simplicity. 
foreach my $sample_id (@sample_table) {
    $sample_count++;
    my $local_time = localtime();
    print "[$local_time] processing sample $sample_id with read-alignment filtering\n";
    my $ref_genome_file = "$base_dir/$ref_genome";
    my $bam_file = "$base_dir/$short_read_mapping_dir/$batch_id/$sample_id/$sample_id.final.bam";
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

    # setup reference genome
    my $min_expected_gc = 0;
    my $max_expected_gc = 1;
    if (-d "$base_dir/$output_dir/ref_genome_preprocessing") {
	print "found pre-calculated ref_geome_preprocessing directory!\n";
	my $min_expected_gc_file = "$base_dir/$output_dir/ref_genome_preprocessing/min_expected_gc.txt";
	my $min_expected_gc_fh = read_file($min_expected_gc_file);
	$min_expected_gc = parse_expected_gc($min_expected_gc_fh);
	close($min_expected_gc_fh);
	my $max_expected_gc_file = "$base_dir/$output_dir/ref_genome_preprocessing/max_expected_gc.txt";
	my $max_expected_gc_fh = read_file($max_expected_gc_file);
	$max_expected_gc = parse_expected_gc($max_expected_gc_fh);
	close($max_expected_gc_fh);

	print "picking up min_expected_gc and max_expected_gc values\n";
	print "min_expected_gc=$min_expected_gc, max_expected_gc=$max_expected_gc\n";
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
	my $for_FREEC_refseq = "ref.genome.fa";
	my $for_FREEC_refseq_fh = read_file($for_FREEC_refseq);
	my %for_FREEC_refseq = ();
	my @for_FREEC_refseq = ();
	parse_fasta_file($for_FREEC_refseq_fh, \%for_FREEC_refseq, \@for_FREEC_refseq);
	
	my $FREEC_refseq = "FREEC.refseq.fa";
	my $FREEC_refseq_fh = write_file($FREEC_refseq);
	my %FREEC_refseq_chr_dict = ();
	my @FREEC_chr = ();
	# setup chr, chrLength   
	system("mkdir FREEC_refseq_chr");
	my $chr_index = 0;
	foreach my $chr (@for_FREEC_refseq) {
	    $chr_index++;
	    my $FREEC_chr = "chr". $chr_index;
	    $FREEC_refseq_chr_dict{'FREEC_to_original'}{$chr_index} = $chr;
	    $FREEC_refseq_chr_dict{'original_to_FREEC'}{$chr} = $FREEC_chr;
	    push @FREEC_chr, $FREEC_chr;
	    print $FREEC_refseq_fh ">$FREEC_chr\n$for_FREEC_refseq{$chr}\n";
	    my $FREEC_chr_fa = "./FREEC_refseq_chr/${FREEC_chr}.fa";
	    my $FREEC_chr_fa_fh = write_file($FREEC_chr_fa);
	    print $FREEC_chr_fa_fh ">$FREEC_chr\n$for_FREEC_refseq{$chr}\n";
	    close $FREEC_chr_fa_fh;
	}
	system("$samtools_dir/samtools faidx FREEC.refseq.fa");
	# determine GC range for FREEC
	my $lower_quantile = 15;
	my $upper_quantile = 100 - $lower_quantile;
	system("$bedtools_dir/bedtools makewindows -g FREEC.refseq.fa.fai -w $window_size > FREEC.refseq.window.$window_size.bed");
	system("$bedtools_dir/bedtools nuc -fi FREEC.refseq.fa -bed FREEC.refseq.window.$window_size.bed > FREEC.refseq.GC_content.txt");
	my $GC = "FREEC.refseq.GC_content.txt";
	my $GC_fh = read_file($GC);
	my %GC = parse_GC_file($GC_fh, $window_size);
	
	foreach my $chr (@FREEC_chr) {
	    my @GC_chr = ();
	    foreach my $i (sort {$a <=> $b} keys %{$GC{$chr}}) {
		if ($GC{$chr}{$i}{'effective_ratio'} > $min_mappability) {
		    push @GC_chr, $GC{$chr}{$i}{'GC'};
		}
	    }
	    if ((scalar @GC_chr) >= 10) {
		my $gc_stat_chr = Statistics::Descriptive::Full->new();
		$gc_stat_chr->add_data(@GC_chr);
		$gc_stat_chr->sort_data();
		$gc_stat_chr->presorted(1);
		# my $gc_mean_chr = sprintf("%.3f", $gc_stat_chr->mean());
		# my $gc_median_chr = sprintf("%.3f", $gc_stat_chr->median());
		my $gc_min_chr = sprintf("%.3f", $gc_stat_chr->min());
		my $gc_max_chr = sprintf("%.3f", $gc_stat_chr->max());
		# my $gc_stdev_chr = sprintf("%.3f", $gc_stat_chr->standard_deviation());
		my $gc_quantile_lower_chr = $gc_stat_chr->percentile($lower_quantile);
		$gc_quantile_lower_chr = sprintf("%.3f", $gc_quantile_lower_chr);
		my $gc_quantile_upper_chr = $gc_stat_chr->percentile($upper_quantile);
		$gc_quantile_upper_chr = sprintf("%.3f", $gc_quantile_upper_chr);
		if ($gc_quantile_lower_chr > $min_expected_gc) {
		    $min_expected_gc = $gc_quantile_lower_chr;
		}
		if (($gc_quantile_upper_chr < $max_expected_gc) and ($gc_quantile_upper_chr > $min_expected_gc)) {
		    $max_expected_gc = $gc_quantile_upper_chr;
		}
		print "chr=$chr\n";
		print "gc_quantile_upper_chr=$gc_quantile_upper_chr, gc_quantile_lower_chr=$gc_quantile_lower_chr\n";
		print "max_expected_gc=$max_expected_gc, min_expected_gc=$min_expected_gc\n";
		print "\n";
	    }
	}
	
	# write down max_expected_gc and min_expected_gc
	system("echo $max_expected_gc > max_expected_gc.txt");
	system("echo $min_expected_gc > min_expected_gc.txt");
	# mappability calculation by gemtools
	system("$gemtools_dir/gemtools index -t $threads  -i FREEC.refseq.fa -o FREEC.refseq.gem");
	system("$gemtools_dir/gem-mappability -T $threads -I FREEC.refseq.gem -l $raw_read_length -m 0.02 -e 0.02 -o FREEC.refseq");
    }
    my $sample_output_dir = "$base_dir/$output_dir/$sample_id";
    system("mkdir -p $sample_output_dir");
    chdir("$sample_output_dir") or die "cannot change directory to: $!\n";
    system("mkdir tmp");
    ## filter bam file by mapping quality
    system("$samtools_dir/samtools view -b -q $min_mapping_quality $base_dir/$short_read_mapping_dir/$batch_id/$sample_id/$sample_id.final.bam  >$sample_id.filtered.bam");
    # index the filtered.bam file
    system("$samtools_dir/samtools index $sample_id.filtered.bam");

    $local_time = localtime();
    print "[$local_time] processing sample $sample_id with CNV calling\n";

    # scan for aneuploidy with FREEC
    if (-s "$base_dir/$excluded_chr_list") {
    	system("perl $VARATHON_HOME/scripts/run_FREEC_wrapper_lite.pl -r ./../ref_genome_preprocessing/ref.genome.fa -bam $sample_id.filtered.bam -prefix $sample_id -ploidy $ploidy -bedtools $bedtools_dir/bedtools -samtools $samtools_dir/samtools -freec $freec_dir/freec -window $window_size -step $step_size -read_length_for_mappability $raw_read_length -min_mappability $min_mappability -mate_orientation 0 -excluded_chr_list $base_dir/$excluded_chr_list -max_expected_gc $max_expected_gc -min_expected_gc $min_expected_gc -threads $threads");
    } else {
    	system("perl $VARATHON_HOME/scripts/run_FREEC_wrapper_lite.pl -r ./../ref.genome.fa -bam $sample_id.final.bam -prefix $sample_id -ploidy $ploidy -bedtools $bedtools_dir/bedtools -samtools $samtools_dir/samtools -freec $freec_dir/freec -window $window_size -step $step_size -read_length_for_mappability $raw_read_length -min_mappability $min_mappability -mate_orientation 0 -max_expected_gc $max_expected_gc -min_expected_gc $min_expected_gc  -threads $threads");
    }
    system("Rscript --vanilla --slave $VARATHON_HOME/scripts/CNV_segmentation_by_DNAcopy.R --input $sample_id.FREEC.bam_ratio.txt --prefix $sample_id --window $window_size --ploidy $ploidy --genome_fai ./../ref_genome_preprocessing/ref.genome.fa.fai");
    system("perl $VARATHON_HOME/scripts/adjust_FREEC_copynumber_by_DNAcopy_copynumber.pl -i $sample_id.FREEC.bam_ratio.sorted.txt -a $sample_id.FREEC.bam_ratio.sorted.resegmented.lite.txt -o $sample_id.FREEC.bam_ratio.sorted.adjusted.txt");
    if (-s "$sample_id.FREEC.bam_ratio.sorted.adjusted.txt") {
    	$local_time = localtime();
    	print "[$local_time] processing sample $sample_id with CNV calls plotting\n";
    	system("Rscript --vanilla --slave $VARATHON_HOME/scripts/plot_CNV_for_FREEC.R --ploidy $ploidy --genome_fai ./../ref_genome_preprocessing/ref.genome.fa.fai --input $sample_id.FREEC.bam_ratio.sorted.adjusted.txt --output $sample_id.CNV_plot.pdf");
    	system("rm Rplots.pdf");
	$local_time = localtime();
	print "[$local_time] processing sample $sample_id with CNV significance test calculation\n";
    	if (-s "$sample_id.FREEC.bam_ratio.sorted.resegmented.CNVs.txt") {
    	    system("Rscript --vanilla --slave $VARATHON_HOME/scripts/assess_CNV_significance_for_FREEC.R --cnv $sample_id.FREEC.bam_ratio.sorted.resegmented.CNVs.txt --ratio $sample_id.FREEC.bam_ratio.sorted.adjusted.txt --genome_fai ./../ref_genome_preprocessing/ref.genome.fa.fai --output $sample_id.CNV_significance_test.txt");
    	    my $sample_cnv_fh = read_file("$sample_id.CNV_significance_test.txt");
    	    parse_sample_cnv_file($sample_cnv_fh, $sample_id, "ref.refseq.fa", $all_samples_cnv_fh);
    	    close $sample_cnv_fh;
    	} else {
    	    system("echo -e \"chr\tstart\tend\tcopy_number\tstatus\tMWU_test_p_value\tKS_test_p_value\" > $sample_id.CNV_significance_test.txt");
    	}
    } else {
    	system(" echo \"Exception encountered for FREEC! Exit! ...\" > $sample_id.final.bam.no_FREEC.txt");
    }
    system("rm -r tmp");
    $local_time = localtime();
    print "[$local_time] finishing short-read-based CNV calling for sample $sample_id \n";
  
    # remove old files
    if ($debug eq "no") {
	# system("rm ref.genome.fa");
	system("rm $sample_id.filtered.bam");
	system("rm $sample_id.filtered.bam.bai");
	system("rm for_CNV.bam");
	system("rm FREEC.bam");
	system("rm FREEC.bam.header.old.sam");
	system("rm FREEC.bam.header.new.sam");
	# system("rm FREEC.refseq.log");
	# system("rm for_CNV.refseq.fa");
	# system("rm for_CNV.refseq.fa.fai");
	# system("rm FREEC.refseq.fa");
	# system("rm FREEC.refseq.fa.fai");
	# system("rm -r FREEC_refseq_chr");
	# system("rm FREEC.refseq.gem");
	# system("rm FREEC.refseq.mappability");
    }
    chdir("./../") or die "cannot change directory to: $!\n";
}

close $all_samples_cnv_fh;

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

sub parse_fasta_file {
    my ($fh, $input_hashref, $input_arrayref) = @_;
    my $seq_name = "";
    while (<$fh>) {
	chomp;
	if (/^\s*$/) {
	    next;
	} elsif (/^\s*#/) {
	    next;
	} elsif (/^>(.*)/) {
    $seq_name = $1;
    push @$input_arrayref, $seq_name;
    $$input_hashref{$seq_name} = "";
	} else {
	    $$input_hashref{$seq_name} .= $_;
	}
    }
}

sub parse_GC_file {
    my ($fh, $window_size) = @_;
    my %GC= ();
    my $i = 0;
    while (<$fh>) {
        chomp;
        /^#/ and next;
	/^\s*$/ and next;
        my ($chr, $window_start, $window_end, $AT_pct, $GC_pct, $A_count, $C_count, $G_count, $T_count, $N_count, $other_count, $window_seq_length) = split /\s+/, $_;
	if (not exists $GC{$chr}) {
	    $i = 0;
	} else {
	    $i++;
	}

	my $AT_count = $A_count + $T_count;
	my $GC_count = $G_count + $C_count;
	$GC{$chr}{$i}{'start'} = $window_start;
	$GC{$chr}{$i}{'end'} = $window_end;
	$GC{$chr}{$i}{'effective_window_size'} = $AT_count + $GC_count;
	$GC{$chr}{$i}{'effective_ratio'} = sprintf("%.3f", $GC{$chr}{$i}{'effective_window_size'}/$window_size);
	if ($GC{$chr}{$i}{'effective_window_size'} > 0) {
	    $GC{$chr}{$i}{'GC'} = $GC_count/$GC{$chr}{$i}{'effective_window_size'};
	} else {
	    $GC{$chr}{$i}{'GC'} = "NA";
	}
    }
    return %GC;
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


sub parse_expected_gc {
    my $fh = shift @_;
    my $expected_gc;
    while (<$fh>) {
	chomp;
	/^#/ and next;
	/^\s*$/ and next;
	$expected_gc = $_;
	last;
    }
    return $expected_gc;
}

 
