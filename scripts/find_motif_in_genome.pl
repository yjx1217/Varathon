#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;

##############################################################
#  script: find_motif_in_genome.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2018.02.28
#  description: find motifs in a multi-fasta file (e.g. a genome assembly)
#  example: perl find_motif_in_genome.pl -i genome.fa(.gz) -m motif.list -p prefix
##############################################################

my ($input, $motif, $prefix);

GetOptions('input|i:s' => \$input,
	   'motif|m:s' => \$motif,
	   'prefix|p:s' => \$prefix);


my $input_fh = read_file($input);
my %input = ();
my @input = ();
parse_fasta_file($input_fh, \%input, \@input);

my $motif_fh = read_file($motif);
my %motif = parse_list_file($motif_fh); 

my $detail = "$prefix.masking_details.txt";
my $detail_fh = write_file($detail);
my $summary = "$prefix.masking_summary.txt";
my $summary_fh = write_file($summary);

my %match;

foreach my $m (sort keys %motif) {
    print "m=$m\n";
    my $i = 0;
    foreach my $s (@input) {
	$input{$s} = uc $input{$s};
	my $total_match_length = 0;
	while ($input{$s} =~ /$m/g) {
	    $i++;
	    $match{$m}{$i}{'case'} = $&; 
	    # my $start = pos($input{$s}) + 1;
	    # my $end = $start + (length $match{$m}{$i}{'case'}) - 1;
	    # $match{$m}{$i}{'position'} = "$s:$start-$end";
	    my $start = $-[0] + 1;
	    my $end = $+[0];
	    $match{$m}{$i}{'position'} = "$s:$start-$end";
	    $match{$m}{$i}{'length'} = $end - $start + 1;
	}
    }
}

my %stat;
print $detail_fh "motif\tmatch_id\tref\tmatch_start\tmatch_end\tmatch_case\n";
foreach my $m (sort keys %motif) {
    foreach my $i (sort {$a<=>$b} keys %{$match{$m}}) {
	my ($s, $start, $end) = ($match{$m}{$i}{'position'} =~ /(\S+):(\d+)-(\d+)/);
	print $detail_fh "$m\t$i\t$s\t$start\t$end\t$match{$m}{$i}{'case'}\n";
	if (exists $stat{$m}{$s}) {
	    $stat{$m}{$s}{'count'}++;
	    $stat{$m}{$s}{'length'} += ($end - $start +1);

	} else {
	    $stat{$m}{$s}{'count'} = 1;
	    $stat{$m}{$s}{'length'} = $end - $start +1;
	}
    }
}

print $summary_fh "motif\tref\tref_length\tmotif_count\tmatch_length\tmatch_cov\n";
foreach my $m (sort keys %motif) {
    foreach my $s (@input) {
	if (exists $stat{$m}{$s}) {
	    my $N_count = $input{$s} =~ tr/N/N/;
	    my $length = (length $input{$s}) - $N_count;
	    $stat{$m}{$s}{'cov'} = $stat{$m}{$s}{'length'}/$length;
	    print $summary_fh "$m\t$s\t$length\t$stat{$m}{$s}{'count'}\t$stat{$m}{$s}{'length'}\t$stat{$m}{$s}{'cov'}\n";
	}
    }
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
    if ($file =~ /.gz$/) {
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
        if (/^\s*$/) {
            next;
	} elsif (/^#/) {
            next;
	} else {
	    my $line = $_;
	    if (exists $list{$line}) {
		$list{$line}++;
	    } else {
		$list{$line} = 1;
	    }
	}
    }
    return %list;
}


sub revcom {
    my $seq = shift @_;
    my $seq_revcom = reverse $seq;
    $seq_revcom =~ tr/ATGCNatgcn/TACGNtacgn/;
    return $seq_revcom;
}
