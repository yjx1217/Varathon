#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Statistics::Descriptive;

##############################################################
#  script: adjust_FREEC_copynumber_by_DNAcopy_copynumber.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2018.09.12
#  description: adjust FREEC copy number by the copy number assigned by DNAcopy
#  example: adjust_FREEC_copy_number_by_DNAcopy_resegmentation.pl -i FREEC_ratio.txt -o FREEC_ratio.adjusted.txt -a DNAcopy.resegmentation.txt
##############################################################

my ($input, $output, $adjusted_copy_number);

GetOptions('input|i:s' => \$input,
           'output|o:s' => \$output,
	   'adjust|a:s' => \$adjusted_copy_number);

my $input_fh = read_file($input);
my $output_fh = write_file($output);
my $adjusted_copy_number_fh = read_file($adjusted_copy_number);
my %adjusted_copy_number = parse_adjusted_copy_number_file($adjusted_copy_number_fh);
close $adjusted_copy_number_fh;


while (<$input_fh>) {
    chomp;
    /^#/ and next;
    /^\s*$/ and next;
    if (/^Chromosome\tStart/) {
	print $output_fh "$_\n";
    } else {
	my ($chr, $start, $end, $ratio, $median_ratio, $copy_number) = split /\t/, $_;
	my $flag = 0;
	if ($ratio ne "-1") {
	    foreach my $i (sort {$a <=> $b} keys %adjusted_copy_number) {
		if ($chr eq $adjusted_copy_number{$i}{'chr'}) {
		    if (($start >= $adjusted_copy_number{$i}{'start'}) and ($end <= $adjusted_copy_number{$i}{'end'})) {
			$copy_number = $adjusted_copy_number{$i}{'copy_number'};
			$flag = 1;
			last;
		    }
		}
	    }
	} 
	if ($flag == 0) {
	    $copy_number = "NA";
	}
	print $output_fh "$chr\t$start\t$end\t$ratio\t$median_ratio\t$copy_number\n";
    }
}

close $input_fh;
close $output_fh;
    



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
    if ($file =~ /.gz$/) {
        open($fh, "| gzip -c >$file") or die "can't open $file\n";
    } else {
        open($fh, ">$file") or die "can't open $file\n";
    }
    return $fh;
}  

sub parse_adjusted_copy_number_file {
    my $fh = shift @_;
    my %copy_number = ();
    my $index = 0;
    while (<$fh>) {
	chomp;
	/^#/ and next;
	/^\s*$/ and next;
	/^Chromosome\tStart/ and next;
	my ($chr, $start, $end, $ratio, $copy_number) = split /\t/, $_;
	$index++;
	$copy_number{$index}{'chr'} = $chr;
	$copy_number{$index}{'start'} = $start;
	$copy_number{$index}{'end'} = $end;
	$copy_number{$index}{'ratio'} = $ratio;
	$copy_number{$index}{'copy_number'} = $copy_number;
    }
    return %copy_number;
}
