#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: coverage2read_number.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2019.04.14
#  description: calculate the number of reads needed to coverage the input genome with specified coverage
#  example: perl coverage2read_number.pl -i input_genome.fa(.gz) -coverage 50 -mean_read_length 4500 -o output_prefix.read_number.txt
##############################################################

my ($input, $output, $coverage, $mean_read_length);
$output = "coverage2read_number.out.txt";
GetOptions('input|i:s' => \$input, # input genome fasta file
	   'coverage|c:f' => \$coverage,
	   'mean_read_length|l:i' => \$mean_read_length,
           'output|o:s' => \$output);

my $input_fh = read_file($input);
my @input = ();
my %input = ();
parse_fasta_file($input_fh, \%input, \@input);

my $genome_size = 0;
foreach my $chr (@input) {
    $genome_size += (length $input{$chr});
}

my $read_number = int ($genome_size * $coverage / $mean_read_length + 0.5);
my $output_fh = write_file($output);
print $output_fh "$read_number\n";

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

sub parse_fasta_file {
    my ($fh, $seq_hashref, $seq_arraryref) = @_;
    my $seq_id = "";
    while (<$fh>) {
        chomp;
        if (/^\s*$/) {
            next;
        } elsif (/^\s*\#/) {
            next;
        } elsif (/^>(\S+)/) {
            $seq_id = $1;
            $$seq_hashref{$seq_id} = "";
            push @$seq_arraryref, $seq_id;
        } else {
            $$seq_hashref{$seq_id} .= $_;
        }
    }
}
