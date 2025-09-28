#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Statistics::Descriptive;

##############################################################
#  script: run_FREEC_wrapper_lite.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2023.05.03
#  description: enable more automatic processing of FREEC (without reference preprocessing)
#  example: 
#   perl run_FREEC_wrapper.pl \
#   -refseq ref.refseq.fa(.gz) \
#   -ploidy 2 \
#   -threads 4 \
#   -samtools <path to the samtools binary> \
#   -bedtools <path_to the bedtools_binary> \
#   -freec <path to the freec binary> \
#   -read_length_for_mappability 100 \
#   -window 10000 \
#   -step 2500 \
#   -mate_orientation 0 \
#   -bam input.bam \
#   -excluded_chr_list excluded_chr_list.txt \
#   -output_prefix output_prefix
##############################################################

my ($refseq, $bam, $prefix, $ploidy, $threads, $bedtools, $samtools, $freec, $window_size, $step_size, $mate_orientation, $read_length_for_mappability, $min_mappability, $min_expected_gc, $max_expected_gc, $excluded_chr_list);

$ploidy = 2;
$threads = 1;
$window_size = 250;
$step_size = $window_size;
$min_mappability = 0.85;
$excluded_chr_list = "";
$min_expected_gc = 0.3; # this value will be automatically adjusted based on the input genome, so no need to change
$max_expected_gc = 0.6; # this value will be automatically adjusted based on the input genome, so no need to change
$mate_orientation = 0; # '0' for sorted bam, 'FR' for unsorted Illumina paired-ends, 'RF' for Illumina mate-pairs
my $telocentromeric = 0;

GetOptions('refseq|r:s' => \$refseq,
	   'bam|bam:s' => \$bam,
	   'prefix|prefix:s' => \$prefix,
	   'ploidy|ploidy:i' => \$ploidy,
	   'threads|t:i' => \$threads,
	   'bedtools|bedtools:s' => \$bedtools,
	   'samtools|samtools:s' => \$samtools,
	   'freec|f:s' => \$freec,
	   'min_mappability|min_map:f' => \$min_mappability,
	   'min_expected_gc|min_gc:f' => \$min_expected_gc,
	   'max_expected_gc|max_gc:f' => \$max_expected_gc,
	   'window|w:i' => \$window_size,
	   'step|s:i' => \$step_size,
	   'read_length_for_mappability|read_length:i' => \$read_length_for_mappability,
	   'mate_orientation|mo:s' => \$mate_orientation,
	   'excluded_chr_list|e:s' => \$excluded_chr_list);

my $refseq_fh = read_file($refseq);
my %refseq = ();
my @refseq = ();
parse_fasta_file($refseq_fh, \%refseq, \@refseq);

my @excluded_chr_list_regexp = "";
if ($excluded_chr_list ne "") {
    print "excluded_chr_list = $excluded_chr_list\n";
    my $excluded_chr_list_fh = read_file($excluded_chr_list);
    my %excluded_chr = parse_list_file($excluded_chr_list_fh);
    foreach my $chr (sort keys %excluded_chr) {
	push @excluded_chr_list_regexp, "/$chr/{d;b}";
    }

    my @refseq_filtered = ();
    foreach my $chr (@refseq) {
	if (not exists $excluded_chr{$chr}) {
	    push @refseq_filtered, $chr;
	} else {
	    delete $refseq{$chr};
	}
    }
    @refseq = @refseq_filtered;
}

# my $FREEC_refseq = "FREEC.refseq.fa";
# my $FREEC_refseq_fh = write_file($FREEC_refseq);
# my @FREEC_refseq = ();
my %FREEC_refseq_chr_dict = ();

# # setup chr, chrLength
# system("mkdir FREEC_refseq_chr");

my $chr_index = 0;
foreach my $chr (@refseq) {
    $chr_index++;
    my $FREEC_chr = "chr". $chr_index;
    $FREEC_refseq_chr_dict{'FREEC_to_original'}{$chr_index} = $chr;
    $FREEC_refseq_chr_dict{'original_to_FREEC'}{$chr} = $FREEC_chr;
    # print $FREEC_refseq_fh ">$FREEC_chr\n$refseq{$chr}\n";
    # my $FREEC_chr_fa = "./FREEC_refseq_chr/${FREEC_chr}.fa";
    # my $FREEC_chr_fa_fh = write_file($FREEC_chr_fa);
    # print $FREEC_chr_fa_fh ">$FREEC_chr\n$refseq{$chr}\n";
    # close $FREEC_chr_fa_fh;
}
# system("$samtools faidx FREEC.refseq.fa");

# # determine GC range for FREEC

# my %GC = ();
# my $lower_quantile = 15;
# my $upper_quantile = 100 - $lower_quantile;
# foreach my $chr (@refseq) {
#     my $remainder = $refseq{$chr};
#     my $i = 0;
#     my @GC_chr = ();
#     while ((length $remainder) >= $window_size) {
# 	$i++;
# 	my $window_start = $step_size * ($i - 1) + 1;
# 	my $window_end = $window_start + $window_size - 1;
# 	my $subseq = substr $remainder, 0, $window_size;
# 	$subseq = uc $subseq;
# 	$GC{$chr}{$i}{'start'} = $window_start;
# 	$GC{$chr}{$i}{'end'} = $window_end;
# 	my $GC_count  = () = $subseq =~ /(G|C)/g;
# 	my $AT_count  = () = $subseq =~ /(A|T)/g;
# 	$GC{$chr}{$i}{'effective_window_size'} = $AT_count + $GC_count;
# 	$GC{$chr}{$i}{'effective_ratio'} = sprintf("%.3f", $GC{$chr}{$i}{'effective_window_size'}/$window_size);
# 	my $GC_window;
# 	if ($GC{$chr}{$i}{'effective_window_size'} > 0) {
# 	    $GC{$chr}{$i}{'GC'} = $GC_count/$GC{$chr}{$i}{'effective_window_size'};
# 	    $GC_window = sprintf("%.3f", $GC{$chr}{$i}{'GC'});
# 	} else {
# 	    $GC{$chr}{$i}{'GC'} = "NA";
# 	    $GC_window = "NA";
# 	}
# 	if (($AT_count + $GC_count)/$window_size >= $min_mappability) {
# 	    push @GC_chr, $GC{$chr}{$i}{'GC'};
# 	}
# 	$remainder = substr $remainder, $step_size;
#     }
#     if ((scalar @GC_chr) >= 10) {
#     my $gc_stat_chr = Statistics::Descriptive::Full->new();
#     $gc_stat_chr->add_data(@GC_chr);
#     $gc_stat_chr->sort_data();
#     $gc_stat_chr->presorted(1);
#     # my $gc_mean_chr = sprintf("%.3f", $gc_stat_chr->mean());
#     # my $gc_median_chr = sprintf("%.3f", $gc_stat_chr->median());
#     my $gc_min_chr = sprintf("%.3f", $gc_stat_chr->min());
#     my $gc_max_chr = sprintf("%.3f", $gc_stat_chr->max());
#     # my $gc_stdev_chr = sprintf("%.3f", $gc_stat_chr->standard_deviation());
#     my $gc_quantile_lower_chr = $gc_stat_chr->percentile($lower_quantile);
#     $gc_quantile_lower_chr = sprintf("%.3f", $gc_quantile_lower_chr);
#     my $gc_quantile_upper_chr = $gc_stat_chr->percentile($upper_quantile);
#     $gc_quantile_upper_chr = sprintf("%.3f", $gc_quantile_upper_chr);
#     if ($gc_quantile_lower_chr < $min_expected_gc) {
# 	$min_expected_gc = $gc_quantile_lower_chr;
#     }
#     if ($gc_quantile_upper_chr > $max_expected_gc) {
# 	$max_expected_gc = $gc_quantile_upper_chr;
#     }
#   }
# }

print "min_expected_gc=$min_expected_gc, max_expected_gc=$max_expected_gc\n";
print "excluded_chr_list = $excluded_chr_list\n";
# filter bam

my $excluded_chr_list_regexp = join ";", @excluded_chr_list_regexp;
$excluded_chr_list_regexp = "'" . $excluded_chr_list_regexp . "'";
print "excluded_chr_list_regexp = $excluded_chr_list_regexp\n";
system("$samtools view -h $bam |sed $excluded_chr_list_regexp |$samtools view -b -h - > for_CNV.bam");

# reheader bam file for FREEC
my $FREEC_bam_header_old = "FREEC.bam.header.old.sam";
system("$samtools view -H for_CNV.bam > $FREEC_bam_header_old");
my $FREEC_bam_header_old_fh = read_file($FREEC_bam_header_old);
my $FREEC_bam_header_new = "FREEC.bam.header.new.sam";
my $FREEC_bam_header_new_fh = write_file($FREEC_bam_header_new);
reformat_FREEC_bam_header($FREEC_bam_header_old_fh, $FREEC_bam_header_new_fh, \%FREEC_refseq_chr_dict);
close $FREEC_bam_header_old_fh;
close $FREEC_bam_header_new_fh;

system("$samtools reheader FREEC.bam.header.new.sam for_CNV.bam > FREEC.bam");

# generate config file for FREEC
my $FREEC_config = "FREEC.config.txt";
my $FREEC_config_fh = write_file($FREEC_config);

print $FREEC_config_fh "###For more options see: http://boevalab.com/FREEC/tutorial.html#CONFIG ###\n[general]\n";
print $FREEC_config_fh "bedtools = $bedtools\n";
print $FREEC_config_fh "samtools = $samtools\n";
print $FREEC_config_fh "maxThreads = $threads\n";
print $FREEC_config_fh "chrLenFile = ./../ref_genome_preprocessing/FREEC.refseq.fa.fai\n";
print $FREEC_config_fh "chrFiles = ./../ref_genome_preprocessing/FREEC_refseq_chr\n";
print $FREEC_config_fh "telocentromeric = 0\n";
print $FREEC_config_fh "ploidy = $ploidy\n";
print $FREEC_config_fh "gemMappabilityFile = ./../ref_genome_preprocessing/FREEC.refseq.mappability\n";
print $FREEC_config_fh "minMappabilityPerWindow = $min_mappability\n";
print $FREEC_config_fh "minExpectedGC = $min_expected_gc\n";
print $FREEC_config_fh "maxExpectedGC = $max_expected_gc\n";
print $FREEC_config_fh "window = $window_size\n";
print $FREEC_config_fh "step = $step_size\n";
print $FREEC_config_fh "breakPointThreshold = 1.2\n";
print $FREEC_config_fh "breakPointType = 2\n";
print $FREEC_config_fh "[sample]\n";
print $FREEC_config_fh "inputFormat = BAM\n";
print $FREEC_config_fh "mateFile = FREEC.bam\n";
print $FREEC_config_fh "mateOrientation = $mate_orientation\n";
close $FREEC_config_fh;

# run FREEC
system("$freec -conf FREEC.config.txt");

my $FREEC_bam_ratio_old = "FREEC.bam_ratio.txt";
if (-s $FREEC_bam_ratio_old) {
    my $FREEC_bam_ratio_old_fh = read_file($FREEC_bam_ratio_old);
    my $FREEC_bam_ratio_new = "$prefix.FREEC.bam_ratio.txt";
    my $FREEC_bam_ratio_new_fh = write_file($FREEC_bam_ratio_new);
    reformat_FREEC_bam_ratio($FREEC_bam_ratio_old_fh, $FREEC_bam_ratio_new_fh, $step_size, \%FREEC_refseq_chr_dict);
    my $FREEC_bam_CNVs_old = "FREEC.bam_CNVs";
    my $FREEC_bam_CNVs_old_fh = read_file($FREEC_bam_CNVs_old);
    my $FREEC_bam_CNVs_new = "$prefix.FREEC.bam_CNVs.txt";
    my $FREEC_bam_CNVs_new_fh = write_file($FREEC_bam_CNVs_new);
    reformat_FREEC_bam_CNVs($FREEC_bam_CNVs_old_fh, $FREEC_bam_CNVs_new_fh, \%FREEC_refseq_chr_dict);
}


sub read_file {
    my $file = shift @_;
    my $fh;
    if ($file =~ /\.gz$/) {
        open($fh, "gunzip -c $file |") or die "can't open pipe to $file";
    } else {
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

sub parse_list_file {
    my $fh = shift @_;
    my %list = ();
    while (<$fh>) {
	chomp;
	/^#/ and next;
	/^\s*$/ and next;
	if (exists $list{$_}) {
	    $list{$_}++;
	} else {
	    $list{$_} = 1;
	}
    }
    return %list
}

sub reformat_FREEC_bam_header {
    my ($FREEC_bam_header_old_fh, $FREEC_bam_header_new_fh, $FREEC_refseq_chr_dict_hashref) = @_;
    while (<$FREEC_bam_header_old_fh>) {
	chomp;
	/^#/ and next;
	/^\s*$/ and next;
	if ($_ =~ /^\@SQ\s+SN:(\S+)/) {
	    my $old_chr = $1;
	    my $new_chr = $$FREEC_refseq_chr_dict_hashref{'original_to_FREEC'}{$old_chr};
	    $_ =~ s/SN:$old_chr/SN:$new_chr/;
	}
	print $FREEC_bam_header_new_fh "$_\n";
    }
}


sub reformat_FREEC_bam_ratio {
    my ($FREEC_bam_ratio_old_fh, $FREEC_bam_ratio_new_fh, $step_size, $FREEC_refseq_chr_dict_hashref) = @_;
    while (<$FREEC_bam_ratio_old_fh>) {
	chomp;
	/^#/ and next;
	/^\s*$/ and next;
	if (/^Chromosome\tStart/) {
	    print $FREEC_bam_ratio_new_fh "Chromosome\tStart\tEnd\tRatio\tMedianRatio\tCopyNumber\n";
	} else {
	    my ($old_chr, $start, $ratio, $median_ratio, $copy_number) = split /\t/, $_;
	    my $new_chr = $$FREEC_refseq_chr_dict_hashref{'FREEC_to_original'}{$old_chr};
	    my $end = $start + $step_size - 1;
	    print $FREEC_bam_ratio_new_fh "$new_chr\t$start\t$end\t$ratio\t$median_ratio\t$copy_number\n";
	}
    }
}

sub reformat_FREEC_bam_CNVs {
    my ($FREEC_bam_CNVs_old_fh, $FREEC_bam_CNVs_new_fh, $FREEC_refseq_chr_dict_hashref) = @_;
    while (<$FREEC_bam_CNVs_old_fh>) {
	chomp;
	/^#/ and next;
	/^\s*$/ and next;
	my ($old_chr, $start, $end, $copy_number, $annotation) = split /\t/, $_;
	my $new_chr = $$FREEC_refseq_chr_dict_hashref{'FREEC_to_original'}{$old_chr};
	print $FREEC_bam_CNVs_new_fh "$new_chr\t$start\t$end\t$copy_number\t$annotation\n";
    }
}

