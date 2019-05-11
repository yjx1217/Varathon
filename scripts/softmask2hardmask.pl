#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: softmask2hardmask.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2018.01.20
#  description: convert a softmasked fasta file into a hardmasked one
#  example: perl softmask2hardmask.pl -i input.fa(.gz)  -o output.fa(.gz) 
##############################################################

my ($input, $output);

GetOptions('input|i:s' => \$input,
	   'output|o:s' => \$output);

my $input_fh = read_file($input);
my %input = ();
my @input = ();
parse_fasta_file($input_fh, \%input, \@input);
close $input_fh;

my $output_fh = write_file($output);

foreach my $id (@input) {
    $input{$id} =~ tr/[atgcn]/N/;
    print $output_fh ">$id\n$input{$id}\n";
}
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


