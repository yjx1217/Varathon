#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: rescale_quality_score_for_clair_vcf.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2020.05.08
#  description: rescaling the quality score computed by clair
#  example: perl rescale_quality_score_for_clair_vcf.pl -i input.vcf(.gz) -o output.vcf(.gz) (-q 30)
##############################################################

my ($input, $output, $qual);
$qual = 0;
GetOptions('input|i:s' => \$input,
	   'output|o:s' => \$output,
	   'qual|q:i' => \$qual);


my $input_fh = read_file($input);
my $output_fh = write_file($output);

while (<$input_fh>) {
    chomp;
    /^\s*$/ and next;
    if (/^#/) {
	print $output_fh "$_\n";
    } else {
	my @line = split /\t/, $_;
	my $old_quality = $line[5];
	my $new_quality = sqrt($old_quality);
	$line[5] = int($new_quality + 0.5);
	if ($line[5] >= $qual) {
            my $line = join "\t", @line;
            print $output_fh "$line\n";
        }
    }
}


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
