#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: filter_vcf_by_bed.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2020.05.14
#  description: filter vcf (for SNP/INDEL/SV) by regions defined in a bed file
#  example: perl filter_vcf_by_bed.pl -i input.vcf(.gz) -b regions.bed -o output.vcf(.gz)
##############################################################

my ($input, $bed, $output);

GetOptions('input|i:s' => \$input,
	   'bed|b:s' => \$bed,
	   'output|o:s' => \$output);


my $input_fh = read_file($input);
my $bed_fh = read_file($bed);
my %bed = parse_bed_file($bed_fh);

my $input_count = 0;
my $output_count = 0;
my $output_fh = write_file($output);

while (<$input_fh>) {
    chomp;
    (/^\s*$/) and next;
    if (/^#/) {
	print $output_fh "$_\n";
    } else {
	my ($chr, $start, $id, $ref_allele, $alt_allele, $qual, $filter, $info) = split /\t/, $_;
	$input_count++;
	my $end;
	if ($info =~ /END=([^;]+)/) {
	    $end = $1;
	} else {
	    $end = abs(length($alt_allele) - length($ref_allele)) + $start;
	}
	my $flag = 0;
	if (not exists $bed{$chr}) {
	    next;
	} else {
	    foreach my $s (sort {$a <=> $b} keys  %{$bed{$chr}}) {
		foreach my $e (sort {$b <=> $a} keys %{$bed{$chr}{$s}}) {
		    if (($start >= $s) and ($end <= $e)) {
			$flag = 1;
			last;
		    }
		}
		if ($flag == 1) {
		    last;
		}
	    }
	    if ($flag == 1) {
		print $output_fh "$_\n";
		$output_count++;
	    }
	}
    }
}

print "input_count: $input_count\n";
print "output_count: $output_count\n";



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

sub parse_bed_file {
    my $fh = shift @_;
    my %bed = ();
    while (<$fh>) {
	chomp;
	(/^\s*$/) and next;
        (/^#/) and next;
	my @line = split /\t/, $_;
	my $chr = $line[0];
	my $start = $line[1] + 1;
	my $end = $line[2];
	$bed{$chr}{$start}{$end} = "$chr:$start-$end";
    }
    return %bed;
}
