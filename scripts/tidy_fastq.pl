#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: tidy_fastq.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2020.05.13
#  description: strip off other tags in the sequence id line that might cause problems for downstream analysis.
#  example: perl tidy_fastq.pl -i input.fastq(.gz) -o output.fastq(.gz)
##############################################################

my ($input, $output);

GetOptions('input|i:s' => \$input,
	   'output|o:s' => \$output);


my $input_fh = read_file($input);
my $count = 0;
my $output_fh = write_file($output);

while (<$input_fh>) {
    chomp;
    if ($. % 4 == 1) {
	my ($id) = ($_ =~ /^\@(\S+)/);
	print $output_fh "\@$id\n";
	$count++;
    } elsif ($. % 4 == 3) {
	print $output_fh "+\n";
    } else {
	print $output_fh "$_\n";
    }
}

print "\nDone! Processed $count reads in total!\n\n";

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


