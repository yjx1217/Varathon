#!/usr/bin/perl
use warnings FATAL => 'all';
use strict;
use Getopt::Long qw(GetOptions);
use Pod::Usage qw(pod2usage);
use List::Util qw(sum min max);
# use Data::Dumper;

##############################################################
#  script: vcf_benchmarker.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2020.06.21
#  description: compare the reference and query vcf files to determine the precision, recall, and F1 score of the query variant calling result by assuming the reference vcf file as the ground truth.
##############################################################

###################
# input parameters
###################

# General options
my $help;
my $man;
my $version;
# Options for the reference genome
my $refseq;
# my $excluded_chr_list;
my $vcf_for_sv = "no";
my $max_offset = 10;
my $snp_indel_qual_filter = 30;
my $sv_qual_filter = "LowQual;q5";
my $ref_vcf;
my $query_vcf;
# Options for the output file prefix
my $prefix = "output_prefix";

GetOptions(
    'help|h|?' => \$help,
    'man|m' => \$man,
    'version|ver' => \$version,
    'vcf_for_sv:s' => \$vcf_for_sv,
    'ref_vcf|r:s' => \$ref_vcf,
    'query_vcf|q:s' => \$query_vcf,
    'snp_indel_qual_filter:s' => \$snp_indel_qual_filter,
    'sv_qual_filter:s' => \$sv_qual_filter,
    'max_offset|d:f' => \$max_offset,
    'prefix|p:s' => \$prefix,
    );

## Option check and print out Usage if needed.
# If usage (or version number) was explicitly requested, print usage (or version number).

if (defined $help) {
  pod2usage(-verbose => 1);
}

if (defined $man) {
  pod2usage(-verbose => 2);
}

if (defined $version) {
  pod2usage(-verbose => 99, -sections => qw(NAME|DESCRIPTION|VERSION));
}

if (not defined $ref_vcf) {
  pod2usage(-message => "Mandatory argument '-ref_vcf' is missing! Exit!",
  -exitval => 1,
  -verbose => 1);
} elsif (not -e $ref_vcf) {
  print "!!! Error! The defined input file $ref_vcf does not exist!\n";
  print "!!! Exit!\n";
  die;
}

if (not defined $query_vcf) {
  pod2usage(-message => "Mandatory argument '-query_vcf' is missing! Exit!",
  -exitval => 1,
  -verbose => 1);
} elsif (not -e $query_vcf) {
  print "!!! Error! The defined input file $query_vcf does not exist!\n";
  print "!!! Exit!\n";
  die;
}

if ($vcf_for_sv eq "no") {
    print "\n";
    print "The variants described in the input vcf files are SNPs/INDELs.\n";
    print "Only consider perfect matches regarding both the position and the alternative allele.\n";
    print "The default variant quality filter $snp_indel_qual_filter will be used when parsing the input VCF files. Records with undefined quality value '.' will not be affected by this filter. Please use '-snp_indel_qual_filter' to redefine it if needed.\n";
    my $ref_vcf_fh = read_file($ref_vcf);
    print "\nParsing the reference VCF file: $ref_vcf\n";
    my %ref_vcf = parse_simple_vcf_file($ref_vcf_fh, $snp_indel_qual_filter);
    close $ref_vcf_fh;
    my $ref_vcf_variant_count = count_variants_for_simple_vcf(\%ref_vcf);
    my $query_vcf_fh = read_file($query_vcf);
    print "Parsing the query VCF file: $query_vcf\n";
    my %query_vcf = parse_simple_vcf_file($query_vcf_fh, $snp_indel_qual_filter);
    close $query_vcf_fh;
    my $query_vcf_variant_count = count_variants_for_simple_vcf(\%query_vcf);
    my %intersect_vcf = simple_vcf_intersect(\%ref_vcf, \%query_vcf);
    my %ref_only_vcf = simple_vcf_uniq(\%intersect_vcf, \%ref_vcf);
    my %query_only_vcf = simple_vcf_uniq(\%intersect_vcf, \%query_vcf);
    
    my $intersect_vcf_variant_count = count_variants_for_simple_vcf(\%intersect_vcf);
    # my $ref_only_vcf_variant_count = count_variants_for_simple_vcf(\%ref_only_vcf);
    # my $query_only_vcf_variant_count = count_variants_for_simple_vcf(\%query_only_vcf);
    # print "intersect_count: $intersect_vcf_variant_count, ref_only_count: $ref_only_vcf_variant_count, query_only_count:$query_only_vcf_variant_count\n";

    my $tp = $intersect_vcf_variant_count;
    my $fn = $ref_vcf_variant_count - $intersect_vcf_variant_count;
    my $fp = $query_vcf_variant_count - $intersect_vcf_variant_count;

    my $recall = "NA";
    my $precision = "NA";
    my $f1 = "NA";

    if ($tp + $fn != 0) {
	$recall = $tp/($tp + $fn);
    }
    if ($tp + $fp != 0) {
	$precision = $tp/($tp + $fp);
    }
    if (($recall ne "NA") and ($precision ne "NA")) {
	if (($recall != 0) and ($precision != 0)) {
	    $f1 = 2 * ($recall * $precision)/($recall + $precision);
	}
    }

    my $output_SNP_INDEL_summary = "$prefix.vcf_benchmarker.SNP_INDEL.summary.txt";
    my $output_SNP_INDEL_summary_fh = write_file($output_SNP_INDEL_summary);
    my $output_SNP_INDEL_detail_ref_intersection = "$prefix.vcf_benchmarker.SNP_INDEL.detail_ref_intersection.txt";
    my $output_SNP_INDEL_detail_ref_intersection_fh = write_file($output_SNP_INDEL_detail_ref_intersection);
    my $output_SNP_INDEL_detail_query_intersection = "$prefix.vcf_benchmarker.SNP_INDEL.detail_query_intersection.txt";
    my $output_SNP_INDEL_detail_query_intersection_fh = write_file($output_SNP_INDEL_detail_query_intersection);
    my $output_SNP_INDEL_detail_ref_only = "$prefix.vcf_benchmarker.SNP_INDEL.detail_ref_only.txt";
    my $output_SNP_INDEL_detail_ref_only_fh = write_file($output_SNP_INDEL_detail_ref_only);
    my $output_SNP_INDEL_detail_query_only = "$prefix.vcf_benchmarker.SNP_INDEL.detail_query_only.txt";
    my $output_SNP_INDEL_detail_query_only_fh = write_file($output_SNP_INDEL_detail_query_only);

    print "\n\n";
    print "\n#########################\n";
    print "ref vcf: $ref_vcf\n";
    print "query vcf: $query_vcf\n";
    print "variant type: SNP/INDEL\n";
    print "variant count in ref vcf: $ref_vcf_variant_count\n";
    print "variant count in query vcf: $query_vcf_variant_count\n";
    print "true positive (ref query intersection): $tp\n";
    print "false positive (query only): $fp\n";
    print "false negative (ref only): $fn\n";
    print "recall: $recall\n";
    print "precision: $precision\n";
    print "f1_score: $f1\n";
    print "#########################\n\n";

    print $output_SNP_INDEL_summary_fh "ref vcf: $ref_vcf\n";
    print $output_SNP_INDEL_summary_fh "query vcf: $query_vcf\n";
    print $output_SNP_INDEL_summary_fh "variant type: SNP/INDEL\n";
    print $output_SNP_INDEL_summary_fh "variant count in ref vcf: $ref_vcf_variant_count\n";
    print $output_SNP_INDEL_summary_fh "variant count in query vcf: $query_vcf_variant_count\n";
    print $output_SNP_INDEL_summary_fh "true positive (ref query intersection): $tp\n";
    print $output_SNP_INDEL_summary_fh "false positive (query only): $fp\n";
    print $output_SNP_INDEL_summary_fh "false negative (ref only): $fn\n";
    print $output_SNP_INDEL_summary_fh "recall: $recall\n";
    print $output_SNP_INDEL_summary_fh "precision: $precision\n";
    print $output_SNP_INDEL_summary_fh "f1_score: $f1\n";

    print $output_SNP_INDEL_detail_ref_intersection_fh "#ref_vcf_chr\tref_vcf_start\tref_vcf_end\tref_vcf_original_allele\tref_vcf_alternative_allele\n";
    foreach my $chr (sort keys %intersect_vcf) {
	foreach my $start (sort {$a <=> $b} keys %{$intersect_vcf{$chr}}) {
	    my $end = $intersect_vcf{$chr}{$start}{'end'};
	    print $output_SNP_INDEL_detail_ref_intersection_fh "$chr\t$start\t$end\t$intersect_vcf{$chr}{$start}{'ref_allele'}\t$intersect_vcf{$chr}{$start}{'alt_allele'}\n";
	}
    }

    print $output_SNP_INDEL_detail_query_intersection_fh "#query_vcf_chr\tquery_vcf_start\tquery_vcf_end\tquery_vcf_original_allele\tquery_vcf_alternative_allele\n";
    foreach my $chr (sort keys %intersect_vcf) {
	foreach my $start (sort {$a <=> $b} keys %{$intersect_vcf{$chr}}) {
	    my $end = $intersect_vcf{$chr}{$start}{'end'};
	    print $output_SNP_INDEL_detail_query_intersection_fh "$chr\t$start\t$end\t$intersect_vcf{$chr}{$start}{'ref_allele'}\t$intersect_vcf{$chr}{$start}{'alt_allele'}\n";
	}
    }

    print $output_SNP_INDEL_detail_ref_only_fh "#ref_vcf_chr\tref_vcf_start\tref_vcf_end\tref_vcf_original_allele\tref_vcf_alternative_allele\n";
    foreach my $chr (sort keys %ref_only_vcf) {
	foreach my $start (sort {$a <=> $b} keys %{$ref_only_vcf{$chr}}) {
	    my $end = $ref_only_vcf{$chr}{$start}{'end'};
	    print $output_SNP_INDEL_detail_ref_only_fh "$chr\t$start\t$end\t$ref_only_vcf{$chr}{$start}{'ref_allele'}\t$ref_only_vcf{$chr}{$start}{'alt_allele'}\n";
	}
    }

    print $output_SNP_INDEL_detail_query_only_fh "#query_vcf_chr\tquery_vcf_start\tquery_vcf_end\tquery_vcf_original_allele\tquery_vcf_alternative_allele\n";
    foreach my $chr (sort keys %query_only_vcf) {
	foreach my $start (sort {$a <=> $b} keys %{$query_only_vcf{$chr}}) {
	    my $end = $query_only_vcf{$chr}{$start}{'end'};
	    print $output_SNP_INDEL_detail_query_only_fh "$chr\t$start\t$end\t$query_only_vcf{$chr}{$start}{'ref_allele'}\t$query_only_vcf{$chr}{$start}{'alt_allele'}\n";
	}
    }
} else {
    $sv_qual_filter = "LowQual;q5";
    print "\n";
    print "The variants described in the input vcf files are structural variants (SVs) as suggested by the specified option '-vcf_for_sv $vcf_for_sv'.\n";
    print "The default variant quality filter $sv_qual_filter will be used when parsing the input SV VCF files. Please use '-f' to redefine it if needed.\n";
    print "The default maximal allowance of basepair offset $max_offset (bp) will be used to determine the breakend match in the reference and query VCF files. Please use '-d' to redefine it if needed.\n";
    my $ref_vcf_fh = read_file($ref_vcf);
    print "\nParsing the reference VCF file: $ref_vcf\n";
    my %ref_vcf = parse_SV_vcf_file($ref_vcf_fh, $sv_qual_filter, $max_offset);
    close $ref_vcf_fh;
    #print "ref_vcf:\n";
    #print Dumper(\%ref_vcf);
    my %ref_vcf_variant_count = count_variants_for_SV_vcf(\%ref_vcf);
    my $query_vcf_fh = read_file($query_vcf);
    print "\nParsing the query VCF file: $query_vcf\n";
    my %query_vcf = parse_SV_vcf_file($query_vcf_fh, $sv_qual_filter, $max_offset);
    close $query_vcf_fh;
    #print "query_vcf:\n";
    #print Dumper(\%query_vcf);
    my %query_vcf_variant_count = count_variants_for_SV_vcf(\%query_vcf);
    # print "ref_vcf_breakend_group_count: $ref_vcf_variant_count{'breakend_grp'}\n";
    # print "ref_vcf_event_count: $ref_vcf_variant_count{'event'}\n";
    # print "query_vcf_breakend_group_count: $query_vcf_variant_count{'breakend_grp'}\n";
    # print "query_vcf_event_count: $query_vcf_variant_count{'event'}\n";

    # debugdebug
    my %intersect_vcf_by_ref = SV_vcf_intersect(\%ref_vcf, \%query_vcf, $max_offset);
    my %intersect_vcf_by_query = SV_vcf_intersect(\%query_vcf, \%ref_vcf, $max_offset);
    my %ref_only_vcf = SV_vcf_uniq(\%ref_vcf, \%query_vcf, $max_offset);
    my %query_only_vcf = SV_vcf_uniq(\%query_vcf, \%ref_vcf, $max_offset);
    my %intersect_query_event_id = ();
    foreach my $ref_event_id (sort keys %{$ref_vcf{'event'}}) {
	if (not exists $intersect_vcf_by_ref{'event'}{$ref_event_id}) {
	    $ref_only_vcf{'event'}{$ref_event_id} = $ref_vcf{'event'}{$ref_event_id};
	} else {
	    my $match_query_event_id = $intersect_vcf_by_ref{'event'}{$ref_event_id}{'match_query_event_id'};
	    $intersect_query_event_id{$match_query_event_id}{$ref_event_id} = 1;
	}
    }

    foreach my $query_event_id (sort keys %{$query_vcf{'event'}}) {
	if (not exists $intersect_query_event_id{$query_event_id}) {
	    $query_only_vcf{'event'}{$query_event_id} = $query_vcf{'event'}{$query_event_id};
	}
    }
    #############
    # print "########### intersect_vcf_by_ref\n";
    # print Dumper(%intersect_vcf_by_ref);
    # print "########### intersect_vcf_by_query\n";
    # print Dumper(%intersect_vcf_by_query);
    # print "##########\n";
    my %intersect_vcf_by_ref_variant_count = count_variants_for_SV_vcf(\%intersect_vcf_by_ref);
    my %intersect_vcf_by_query_variant_count = count_variants_for_SV_vcf(\%intersect_vcf_by_query);
    # my %ref_only_vcf_variant_count = count_variants_for_SV_vcf(\%ref_only_vcf);
    # my %query_only_vcf_variant_count = count_variants_for_SV_vcf(\%query_only_vcf);

    # comparison for breakends
    print "\nComparison for breakend groups (merging distance = $max_offset) >> \n";
    print "breakend group count in the input ref vcf: $ref_vcf_variant_count{'breakend_grp'}\n";
    print "breakend group count in the input query vcf: $query_vcf_variant_count{'breakend_grp'}\n";
    my $bp_tp = "NA";
    my $bp_fn = "NA";
    my $bp_fp = "NA";
    my $bp_recall = "NA";
    my $bp_precision = "NA";
    my $bp_f1 = "NA";

    if ($ref_vcf_variant_count{'breakend_grp'} == 0) {
	print "Warnings! No breakpiont was found in the input ref vcf file: $ref_vcf! Please double check!\n"; 
	print "Exit!\n";
    } else {
	if ($query_vcf_variant_count{'breakend_grp'} == 0) {
	    print "Warnings! No breakpiont was found in the input query vcf file: $query_vcf! You might want to double check!\n"; 
	    $bp_tp = 0;
	    $bp_fn = $ref_vcf_variant_count{'breakend_grp'};
	    $bp_fp = 0;
	    $bp_recall = 0;
	    $bp_precision = "NA";
	    $bp_f1 = "NA";
	} else {
	    $bp_tp = $intersect_vcf_by_ref_variant_count{'breakend_grp'};
	    $bp_fn = $ref_vcf_variant_count{'breakend_grp'} - $intersect_vcf_by_ref_variant_count{'breakend_grp'};
	    $bp_fp = $query_vcf_variant_count{'breakend_grp'} - $intersect_vcf_by_query_variant_count{'breakend_grp'};
	
	    if ($bp_tp + $bp_fn != 0) {
		$bp_recall = $bp_tp/($bp_tp + $bp_fn);
	    }
	    if ($bp_tp + $bp_fp != 0) {
		$bp_precision = $bp_tp/($bp_tp + $bp_fp);
	    }
	    if (($bp_recall ne "NA") and ($bp_precision ne "NA")) {
		if (($bp_recall != 0) and ($bp_precision != 0)) {
		    $bp_f1 = 2 * ($bp_recall * $bp_precision)/($bp_recall + $bp_precision);
		}
	    }
	}
    }

    print "\n#########################\n";
    print "ref vcf: $ref_vcf\n";
    print "query vcf: $query_vcf\n";
    print "variant type: SV (breakend group)\n";
    print "maximal matching offset cutoff: $max_offset\n";
    print "breakend group count in ref vcf: $ref_vcf_variant_count{'breakend_grp'}\n";
    print "breakend group count in query vcf: $query_vcf_variant_count{'breakend_grp'}\n";
    print "true positive (ref query intersection): $bp_tp\n";
    print "false positive (query only): $bp_fp\n";
    print "false negative (ref only): $bp_fn\n";
    print "recall: $bp_recall\n";
    print "precision: $bp_precision\n";
    print "f1_score: $bp_f1\n";
    print "########################\n\n";

    my $output_SV_breakend_group_summary = "$prefix.vcf_benchmarker.SV_breakend_group.summary.txt";
    my $output_SV_breakend_group_summary_fh = write_file($output_SV_breakend_group_summary);
    my $output_SV_breakend_group_detail_ref_intersection = "$prefix.vcf_benchmarker.SV_breakend_group.detail_ref_intersection.txt";
    my $output_SV_breakend_group_detail_ref_intersection_fh = write_file($output_SV_breakend_group_detail_ref_intersection);
    my $output_SV_breakend_group_detail_query_intersection = "$prefix.vcf_benchmarker.SV_breakend_group.detail_query_intersection.txt";
    my $output_SV_breakend_group_detail_query_intersection_fh = write_file($output_SV_breakend_group_detail_query_intersection);
    my $output_SV_breakend_group_detail_ref_only = "$prefix.vcf_benchmarker.SV_breakend_group.detail_ref_only.txt";
    my $output_SV_breakend_group_detail_ref_only_fh = write_file($output_SV_breakend_group_detail_ref_only);
    my $output_SV_breakend_group_detail_query_only = "$prefix.vcf_benchmarker.SV_breakend_group.detail_query_only.txt";
    my $output_SV_breakend_group_detail_query_only_fh = write_file($output_SV_breakend_group_detail_query_only);

    print $output_SV_breakend_group_summary_fh "ref vcf: $ref_vcf\n";
    print $output_SV_breakend_group_summary_fh "query vcf: $query_vcf\n";
    print $output_SV_breakend_group_summary_fh "variant type: SV (breakend_group)\n";
    print $output_SV_breakend_group_summary_fh "maximal matching offset cutoff: $max_offset\n";
    print $output_SV_breakend_group_summary_fh "breakend group count in ref vcf: $ref_vcf_variant_count{'breakend_grp'}\n";
    print $output_SV_breakend_group_summary_fh "breakend group count in query vcf: $query_vcf_variant_count{'breakend_grp'}\n";
    print $output_SV_breakend_group_summary_fh "true positive (ref query intersection): $bp_tp\n";
    print $output_SV_breakend_group_summary_fh "false positive (query only): $bp_fp\n";
    print $output_SV_breakend_group_summary_fh "false negative (ref only): $bp_fn\n";
    print $output_SV_breakend_group_summary_fh "recall: $bp_recall\n";
    print $output_SV_breakend_group_summary_fh "precision: $bp_precision\n";
    print $output_SV_breakend_group_summary_fh "f1_score: $bp_f1\n";

    print $output_SV_breakend_group_detail_ref_intersection_fh "#ref_vcf_breakend_group\tref_vcf_breakend_chr\tref_vcf_breakend_start\n";
    foreach my $ref_breakend_grp_id (sort keys %{$intersect_vcf_by_ref{'breakend_grp'}}) {
	foreach my $ref_chr (sort keys %{$ref_vcf{'breakend_grp'}{$ref_breakend_grp_id}{'breakend_hash'}}) {
	    foreach my $ref_start (sort {$a <=> $b} keys %{$ref_vcf{'breakend_grp'}{$ref_breakend_grp_id}{'breakend_hash'}{$ref_chr}}) {
		print $output_SV_breakend_group_detail_ref_intersection_fh "$ref_breakend_grp_id\t$ref_chr\t$ref_start\n";
	    }
	}
	print $output_SV_breakend_group_detail_ref_intersection_fh "\n";
    }

    print $output_SV_breakend_group_detail_query_intersection_fh "#query_vcf_breakend_group\tquery_vcf_breakend_chr\tquery_vcf_breakend_start\n";
    foreach my $ref_breakend_grp_id (sort keys %{$intersect_vcf_by_ref{'breakend_grp'}}) {
	my @match_query_breakend_grp_id = sort keys %{$intersect_vcf_by_ref{'breakend_grp'}{$ref_breakend_grp_id}};
	foreach my $query_breakend_grp_id (@match_query_breakend_grp_id) {
	    foreach my $query_chr (sort keys %{$query_vcf{'breakend_grp'}{$query_breakend_grp_id}{'breakend_hash'}}) {
		foreach my $query_start (sort {$a <=> $b} keys %{$query_vcf{'breakend_grp'}{$query_breakend_grp_id}{'breakend_hash'}{$query_chr}}) {
		    print $output_SV_breakend_group_detail_query_intersection_fh "$query_breakend_grp_id\t$query_chr\t$query_start\n";
		}
	    }
	}
	print $output_SV_breakend_group_detail_query_intersection_fh "\n";
    }

    print $output_SV_breakend_group_detail_ref_only_fh "#ref_vcf_breakend_group\tref_vcf_breakend_chr\tref_vcf_breakend_start\n";
    foreach my $ref_breakend_grp_id (sort keys %{$ref_only_vcf{'breakend_grp'}}) {
	foreach my $ref_chr (sort keys %{$ref_only_vcf{'breakend_grp'}{$ref_breakend_grp_id}{'breakend_hash'}}) {
	    foreach my $ref_start (sort {$a <=> $b} keys %{$ref_only_vcf{'breakend_grp'}{$ref_breakend_grp_id}{'breakend_hash'}{$ref_chr}}) {
		print $output_SV_breakend_group_detail_ref_only_fh "$ref_breakend_grp_id\t$ref_chr\t$ref_start\n";
	    }
	}
	print $output_SV_breakend_group_detail_ref_only_fh "\n";
    }

    print $output_SV_breakend_group_detail_query_only_fh "#query_vcf_breakend_group\tquery_vcf_breakend_chr\tquery_vcf_breakend_start\n";
    foreach my $query_breakend_grp_id (sort keys %{$query_only_vcf{'breakend_grp'}}) {
	foreach my $query_chr (sort keys %{$query_only_vcf{'breakend_grp'}{$query_breakend_grp_id}{'breakend_hash'}}) {
	    foreach my $query_start (sort {$a <=> $b} keys %{$query_only_vcf{'breakend_grp'}{$query_breakend_grp_id}{'breakend_hash'}{$query_chr}}) {
		print $output_SV_breakend_group_detail_query_only_fh "$query_breakend_grp_id\t$query_chr\t$query_start\n";
	    }
	}
	print $output_SV_breakend_group_detail_query_only_fh "\n";
    }

    # comparison for events
    print "\nComparison for events >> \n";
    print "event count in the input ref vcf: $ref_vcf_variant_count{'event'}\n";
    print "event count in the input query vcf: $query_vcf_variant_count{'event'}\n";

    my $event_tp = "NA";
    my $event_fn = "NA";
    my $event_fp = "NA";
    my $event_recall = "NA";
    my $event_precision = "NA";
    my $event_f1 = "NA";

    if ($ref_vcf_variant_count{'event'} == 0) {
	print "Warnings! No event was found in the input ref vcf file: $ref_vcf! Please double check!\n"; 
	print "Exit!\n";
    } else {
	if ($query_vcf_variant_count{'event'} == 0) {
	    print "Warnings! No event was found in the input query vcf file: $query_vcf! You might want to double check!\n"; 
	    $event_tp = 0;
	    $event_fn = $ref_vcf_variant_count{'event'};
	    $event_fp = 0;
	    $event_recall = 0;
	    $event_precision = "NA";
	    $event_f1 = "NA";
	} else {
	    $event_tp = $intersect_vcf_by_ref_variant_count{'event'};
	    $event_fn = $ref_vcf_variant_count{'event'} - $intersect_vcf_by_ref_variant_count{'event'};
	    $event_fp = $query_vcf_variant_count{'event'} - $intersect_vcf_by_query_variant_count{'event'};
	    # print "event_tp=$event_tp, event_fn = $event_fn, event_fp = $event_fp\n";
	    if ($event_tp + $event_fn != 0) {
		$event_recall = $event_tp/($event_tp + $event_fn);
	    }
	    if ($event_tp + $event_fp != 0) {
		$event_precision = $event_tp/($event_tp + $event_fp);
	    }
	    if (($event_recall ne "NA") and ($event_precision ne "NA")) {
		# print "event_recall=$event_recall, event_precision=$event_precision\n";
		if (($event_recall != 0) and ($event_precision != 0)) {
		    $event_f1 = 2 * ($event_recall * $event_precision)/($event_recall + $event_precision);
		}
	    }
	}
    }

    print "\n#########################\n";
    print "ref vcf: $ref_vcf\n";
    print "query vcf: $query_vcf\n";
    print "variant type: SV (event)\n";
    print "maximal matching offset cutoff: $max_offset\n";
    print "event count in ref vcf: $ref_vcf_variant_count{'event'}\n";
    print "event count in query vcf: $query_vcf_variant_count{'event'}\n";
    print "true positive (ref query intersection): $event_tp\n";
    print "false positive (query only): $event_fp\n";
    print "false negative (ref only): $event_fn\n";
    print "recall: $event_recall\n";
    print "precision: $event_precision\n";
    print "f1_score: $event_f1\n";
    print "########################\n\n";


    my $output_SV_event_summary = "$prefix.vcf_benchmarker.SV_event.summary.txt";
    my $output_SV_event_summary_fh = write_file($output_SV_event_summary);
    my $output_SV_event_detail_ref_intersection = "$prefix.vcf_benchmarker.SV_event.detail_ref_intersection.txt";
    my $output_SV_event_detail_ref_intersection_fh = write_file($output_SV_event_detail_ref_intersection);
    my $output_SV_event_detail_query_intersection = "$prefix.vcf_benchmarker.SV_event.detail_query_intersection.txt";
    my $output_SV_event_detail_query_intersection_fh = write_file($output_SV_event_detail_query_intersection);
    my $output_SV_event_detail_ref_only = "$prefix.vcf_benchmarker.SV_event.detail_ref_only.txt";
    my $output_SV_event_detail_ref_only_fh = write_file($output_SV_event_detail_ref_only);
    my $output_SV_event_detail_query_only = "$prefix.vcf_benchmarker.SV_event.detail_query_only.txt";
    my $output_SV_event_detail_query_only_fh = write_file($output_SV_event_detail_query_only);

    print $output_SV_event_summary_fh "ref vcf: $ref_vcf\n";
    print $output_SV_event_summary_fh "query vcf: $query_vcf\n";
    print $output_SV_event_summary_fh "variant type: SV (event)\n";
    print $output_SV_event_summary_fh "maximal matching offset cutoff: $max_offset\n";
    print $output_SV_event_summary_fh "event count in ref vcf: $ref_vcf_variant_count{'event'}\n";
    print $output_SV_event_summary_fh "event count in query vcf: $query_vcf_variant_count{'event'}\n";
    print $output_SV_event_summary_fh "true positive (ref query intersection): $event_tp\n";
    print $output_SV_event_summary_fh "false positive (query only): $event_fp\n";
    print $output_SV_event_summary_fh "false negative (ref only): $event_fn\n";
    print $output_SV_event_summary_fh "recall: $event_recall\n";
    print $output_SV_event_summary_fh "precision: $event_precision\n";
    print $output_SV_event_summary_fh "f1_score: $event_f1\n";

    print $output_SV_event_detail_ref_intersection_fh "#ref_vcf_event_id\tref_vcf_event_name\tref_vcf_event_type\tquery_vcf_event_id\tquery_vcf_event_name\tquery_vcf_event_type\n";
    foreach my $ref_event_id (sort keys %{$intersect_vcf_by_ref{'event'}}) {
	my $ref_event_type = $ref_vcf{'event'}{$ref_event_id}{'sv_event_type'};
	my $ref_event_name = $ref_vcf{'event'}{$ref_event_id}{'sv_event_name'};
	my $query_event_id = $intersect_vcf_by_ref{'event'}{$ref_event_id}{'match_query_event_id'};
	my $query_event_type = $query_vcf{'event'}{$query_event_id}{'sv_event_type'};
	my $query_event_name = $query_vcf{'event'}{$query_event_id}{'sv_event_name'};
	print $output_SV_event_detail_ref_intersection_fh "$ref_event_id\t$ref_event_name\t$ref_event_type\t$query_event_id\t$query_event_name\t$query_event_type\n";
    }

    print $output_SV_event_detail_query_intersection_fh "#query_vcf_event_id\tquery_vcf_event_name\tquery_vcf_event_type\tref_vcf_event_id\tref_vcf_event_name\tref_vcf_event_type\n";
    foreach my $ref_event_id (sort keys %{$intersect_vcf_by_ref{'event'}}) {
	my $ref_event_type = $ref_vcf{'event'}{$ref_event_id}{'sv_event_type'};
	my $ref_event_name = $ref_vcf{'event'}{$ref_event_id}{'sv_event_name'};
	my $query_event_id = $intersect_vcf_by_ref{'event'}{$ref_event_id}{'match_query_event_id'};
	my $query_event_type = $query_vcf{'event'}{$query_event_id}{'sv_event_type'};
	my $query_event_name = $query_vcf{'event'}{$query_event_id}{'sv_event_name'};
	print $output_SV_event_detail_query_intersection_fh "$query_event_id\t$query_event_name\t$query_event_type\t$ref_event_id\t$ref_event_name\t$ref_event_type\n";
    }

    print $output_SV_event_detail_ref_only_fh "#ref_vcf_event_id\tref_vcf_event_name\tref_vcf_event_type\n";
    foreach my $ref_event_id (sort keys %{$ref_only_vcf{'event'}}) {
	my $ref_event_type = $ref_vcf{'event'}{$ref_event_id}{'sv_event_type'};
	my $ref_event_name = $ref_vcf{'event'}{$ref_event_id}{'sv_event_name'};
	print $output_SV_event_detail_ref_only_fh "$ref_event_id\t$ref_event_name\t$ref_event_type\n";
    }

    print $output_SV_event_detail_query_only_fh "#query_vcf_event_id\tquery_vcf_event_name\tquery_vcf_event_type\n";
    foreach my $query_event_id (sort keys %{$query_only_vcf{'event'}}) {
	my $query_event_type = $query_vcf{'event'}{$query_event_id}{'sv_event_type'};
	my $query_event_name = $query_vcf{'event'}{$query_event_id}{'sv_event_name'};
	print $output_SV_event_detail_query_only_fh "$query_event_id\t$query_event_name\t$query_event_type\n";
    }

    # print "################\n";
    # print "ref_only_vcf:\n";
    # print Dumper(%ref_only_vcf);
    # print "################\n";
    # print "query_only_vcf:\n";
    # print Dumper(%query_only_vcf);

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
    open($fh, "| gzip -c >$file") or die "can't open pipe to $file\n";
  } else {
    open($fh, ">$file") or die "can't open $file\n";
  }
  return $fh;
}

sub parse_list_file {
  my $fh = shift @_;
  my %list = ();
  while (<$fh>) {
    chomp;
    /^\s*$/ and next;
    /^#/ and next;
    if (exists $list{$_}) {
      $list{$_}++;
    } else {
      $list{$_} = 1;
    }
  }
  return %list;
}

sub parse_simple_vcf_file {
  my ($fh, $qual_cutoff, $query_type) = @_;
  if (not defined $qual_cutoff) {
      $qual_cutoff = 0;
  }
  my %vcf = ();
  while (<$fh>) {
    chomp;
    /^#/ and next;
    /^\s*$/ and next;
    my ($ref_chr, $ref_start, $variant_id, $ref_allele, $alt_allele, $variant_qual, $variant_filter, $variant_info) = split /\t/, $_;
    if (($variant_filter eq ".") or ($variant_filter eq "PASS")) {
	if (($variant_qual eq ".") or ($variant_qual >= $qual_cutoff)) {
	    my $variant_type;
	    my $ref_allele_length = length $ref_allele;
	    my $alt_allele_length = length $alt_allele;
	    my $ref_end = $ref_start + $ref_allele_length - 1;
	    if (($alt_allele !~ /,/) and ($ref_allele_length ne $alt_allele_length)) {
		$variant_type = "INDEL";
	    } else {
		$variant_type = "SNP";
	    }
	    if ((defined $query_type) and ($query_type ne $variant_type)) {
		next;
	    } else {
		$vcf{$ref_chr}{$ref_start}{'ref_chr'} = $ref_chr;
		$vcf{$ref_chr}{$ref_start}{'ref_start'} = $ref_start;
		$vcf{$ref_chr}{$ref_start}{'ref_end'} = $ref_end;
		$vcf{$ref_chr}{$ref_start}{'ref_allele'} = $ref_allele;
		$vcf{$ref_chr}{$ref_start}{'alt_allele'} = $alt_allele;
		$vcf{$ref_chr}{$ref_start}{'variant_type'} = $variant_type;
		$vcf{$ref_chr}{$ref_start}{'variant_id'} = $variant_id;
		$vcf{$ref_chr}{$ref_start}{'variant_qual'} = $variant_qual;
		$vcf{$ref_chr}{$ref_start}{'variant_info'} = $variant_info;
	    }
	}
    }
  }
  return %vcf;
}

sub count_variants_for_simple_vcf {
    my $vcf_hashref = shift @_;
    my $n = 0;
    foreach my $chr (sort keys %$vcf_hashref) {
	my $i = scalar (keys %{$$vcf_hashref{$chr}});
	$n += $i;
    }
    return $n;
}

sub simple_vcf_intersect {
    my ($ref_vcf_hashref, $query_vcf_hashref) = @_;
    my %intersect_vcf = ();
    foreach my $chr (sort keys %$ref_vcf_hashref) {
        foreach  my $start (sort {$a <=> $b} keys %{$$ref_vcf_hashref{$chr}}) {
            if (exists $$query_vcf_hashref{$chr}{$start}) {
		my $check_ref_allele_flag = 1;
		my $check_alt_allele_flag = 1;
		
                if ($$ref_vcf_hashref{$chr}{$start}{'ref_allele'} eq $$query_vcf_hashref{$chr}{$start}{'ref_allele'}) {
		    $check_ref_allele_flag = 0;
		}
		if (($$ref_vcf_hashref{$chr}{$start}{'alt_allele'} !~ /,/) and $$query_vcf_hashref{$chr}{$start}{'alt_allele'} !~ /,/) {
		    if ($$ref_vcf_hashref{$chr}{$start}{'alt_allele'} eq $$query_vcf_hashref{$chr}{$start}{'alt_allele'}) {
			$check_alt_allele_flag = 0;
		    }
		} elsif (($$ref_vcf_hashref{$chr}{$start}{'alt_allele'} =~ /,/) and $$query_vcf_hashref{$chr}{$start}{'alt_allele'} =~ /,/) {
                    my @ref_alt_allele = split /,/, $$ref_vcf_hashref{$chr}{$start}{'alt_allele'};
		    my @ref_alt_allele_sorted = sort @ref_alt_allele;
		    my $ref_alt_allele_sorted = join ",", @ref_alt_allele_sorted;
                    my @query_alt_allele = split /,/, $$query_vcf_hashref{$chr}{$start}{'alt_allele'};
		    my @query_alt_allele_sorted = sort @query_alt_allele;
		    my $query_alt_allele_sorted = join ",", @query_alt_allele_sorted;
		    if ($ref_alt_allele_sorted eq $query_alt_allele_sorted) {
			$check_alt_allele_flag = 0;
		    }
		} else {
		    $check_alt_allele_flag = 1;
		}
		
		if (($check_ref_allele_flag == 0) and ($check_alt_allele_flag == 0)) {
                    $intersect_vcf{$chr}{$start}{'start'} = $$ref_vcf_hashref{$chr}{$start}{'ref_start'};
                    $intersect_vcf{$chr}{$start}{'end'} = $$ref_vcf_hashref{$chr}{$start}{'ref_end'};
                    $intersect_vcf{$chr}{$start}{'ref_allele'} = $$ref_vcf_hashref{$chr}{$start}{'ref_allele'};
                    $intersect_vcf{$chr}{$start}{'alt_allele'} = $$ref_vcf_hashref{$chr}{$start}{'alt_allele'};
                } 
            }
        }
    }
    return %intersect_vcf;
}

sub simple_vcf_uniq {
    my ($ref_vcf_hashref, $query_vcf_hashref) = @_;
    my %uniq_vcf = ();
    foreach my $chr (sort keys %$query_vcf_hashref) {
        foreach  my $start (sort {$a <=> $b} keys %{$$query_vcf_hashref{$chr}}) {
	    my $flag = 0;
            if (not exists $$ref_vcf_hashref{$chr}{$start}) {
		$flag = 1;
	    } else {
		if (($$ref_vcf_hashref{$chr}{$start}{'alt_allele'} !~ /,/) and $$query_vcf_hashref{$chr}{$start}{'alt_allele'} !~ /,/) {
                    if ($$ref_vcf_hashref{$chr}{$start}{'alt_allele'} ne $$query_vcf_hashref{$chr}{$start}{'alt_allele'}) {
                        $flag = 1;
                    }
                } elsif (($$ref_vcf_hashref{$chr}{$start}{'alt_allele'} =~ /,/) and $$query_vcf_hashref{$chr}{$start}{'alt_allele'} =~ /,/) {
                    my @ref_alt_allele = split /,/, $$ref_vcf_hashref{$chr}{$start}{'alt_allele'};
                    my @ref_alt_allele_sorted = sort @ref_alt_allele;
                    my $ref_alt_allele_sorted = join ",", @ref_alt_allele_sorted;
                    my @query_alt_allele = split /,/, $$query_vcf_hashref{$chr}{$start}{'alt_allele'};
                    my @query_alt_allele_sorted = sort @query_alt_allele;
                    my $query_alt_allele_sorted = join ",", @query_alt_allele_sorted;
                    if ($ref_alt_allele_sorted ne $query_alt_allele_sorted) {
			$flag = 1;
                    }
                } else {
                    $flag = 1;
                }
	    }
	    if ($flag == 1) {
		$uniq_vcf{$chr}{$start}{'start'} = $$query_vcf_hashref{$chr}{$start}{'ref_start'};
		$uniq_vcf{$chr}{$start}{'end'} = $$query_vcf_hashref{$chr}{$start}{'ref_end'};
		$uniq_vcf{$chr}{$start}{'ref_allele'} = $$query_vcf_hashref{$chr}{$start}{'ref_allele'};
		$uniq_vcf{$chr}{$start}{'alt_allele'} = $$query_vcf_hashref{$chr}{$start}{'alt_allele'};
	    }
	}
    }
    return %uniq_vcf;
}

sub add_and_resort {
    my ($a_list, $b_list) = @_;
    # print "a_list=$a_list, b_list=$b_list\n";
    my @a_list;
    my @b_list;
    my %c_list = ();
    if ($a_list =~ /;/) {
	@a_list = split /;/, $a_list;
    } else {
	@a_list = ($a_list);
    }
    if ($b_list =~ /;/) {
	@b_list = split /;/, $b_list;
    } else {
	@b_list = ($b_list);
    }
    foreach my $a (@a_list) {
	if (not exists $c_list{$a}) {
	    $c_list{$a} = 1;
	}
    }
    foreach my $b (@b_list) {
	if (not exists $c_list{$b}) {
	    $c_list{$b} = 1;
	}
    }
    my @c_list = sort keys %c_list;
    my $c_list = join ";", @c_list;
    return $c_list;
}


sub parse_SV_vcf_file {
    my ($fh, $sv_qual_filter, $max_offset) = @_;
    my %sv = ();
    my $record_index = 0;
    my $sv_qual_filter_regexp = $sv_qual_filter;
    $sv_qual_filter_regexp =~ s/;/\|/gi;
    # print "sv_qual_filter=$sv_qual_filter, sv_qual_filter_regexp=$sv_qual_filter_regexp\n";

    while (<$fh>) {
        chomp;
        /^\s*$/ and next;
        /^#/ and next;
        my ($ref_chr, $ref_start, $id, $ref_allele, $alt_allele, $variant_qual, $variant_filter, $variant_info) = split /\t/, $_;
	if ($variant_filter =~ /($sv_qual_filter_regexp)/i) {
            print "Warning! Ignore low quality SV record: $_\n";
            next;
	# } elsif (($variant_filter eq "MapQual") or ($variant_filter eq "")) {
	#    $variant_filter = "PASS";
        } elsif (($variant_filter ne ".") and ($variant_filter ne "PASS")) { 
            print "Warning! Ignore filter failed SV record: $_\n";
	    next;
	} else {
	    $record_index++;
        }
        my $sv_event_type = "NA";
	my $sv_event_name = "NA";
        if ($variant_info =~ /SVTYPE=([^;]+)/) {
            $sv_event_type = $1;
        }
	if ($variant_info =~ /EVENT=([^;]+)/) {
            $sv_event_name = $1;
        } 
        if ($sv_event_type !~ /BND/) {
            my $ref_end;
            if ($variant_info =~ /END=([^;]+)/) {
                $ref_end = $1;
	    } else {
                print "!!! Error! The mandatory tag 'END=' has not been specified in the input vcf for the none BND SV record: $_\n";
                print "!!! $_\n";
                print "!!! Exit!\n";
                die;
            }

	    ($ref_start, $ref_end) = sort {$a <=> $b} ($ref_start, $ref_end);
	    my $ref_start_adj;
	    my $ref_end_adj;
	    if ($sv_event_type eq "DEL") {
		$ref_start_adj = $ref_start + 1;
		$ref_end_adj = $ref_end  + 1;
	    } elsif ($sv_event_type eq "INV") {
		$ref_start_adj = $ref_start - 1;
		$ref_end_adj = $ref_end  + 1;
	    } else {
		$ref_start_adj = $ref_start - 1;
		$ref_end_adj = $ref_end  + 1;
	    }		
	    my @sv_breakend = sort ("$ref_chr:$ref_start", "$ref_chr:$ref_end", "$ref_chr:$ref_start_adj", "$ref_chr:$ref_end_adj");
	    foreach my $breakend (@sv_breakend) {
		my ($breakend_chr, $breakend_start) = split /:/, $breakend;
		if (not exists $sv{'breakend'}{$breakend_chr}{$breakend_start}{'sv_event_type_hash'}{$sv_event_type}) {
		    $sv{'breakend'}{$breakend_chr}{$breakend_start}{'sv_event_type_hash'}{$sv_event_type} = 1;
		}
		if (not exists $sv{'breakend'}{$breakend_chr}{$breakend_start}{'sv_event_name_hash'}{$sv_event_name}) {
		    $sv{'breakend'}{$breakend_chr}{$breakend_start}{'sv_event_name_hash'}{$sv_event_name} = 1;
		}
		foreach my $mate_breakend (@sv_breakend) {
		    if (not exists $sv{'breakend'}{$breakend_chr}{$breakend_start}{'sv_event_mate_breakend_hash'}{$mate_breakend}) {
			$sv{'breakend'}{$breakend_chr}{$breakend_start}{'sv_event_mate_breakend_hash'}{$mate_breakend} = 1;
		    }
		}
	    }
        } elsif ($sv_event_type eq "BND") {
            # see the secion 5 of VCFv4.1 specification (https://samtools.github.io/hts-specs/VCFv4.1.pdf) for the detailed meaning of s, t, and p used below.
            my $s = $ref_allele;
            my $t;
            my $p;
            my $p_relative_strand; # the relative strand of p relative to its original sequence
            my $p_relative_position; # the relative positon of p relative to t: "before_t" or "after_t"
            # print "alt_allele = $alt_allele\n";
            if ($alt_allele =~ /\[$/) {
                $p_relative_strand = "+";
                $p_relative_position = "after_ref_allele";
                ($t, $p) = ($alt_allele =~ /(\S+)\[(\S+)\[/);
            } elsif ($alt_allele =~ /\]$/) {
                $p_relative_strand = "-";
                $p_relative_position = "after_ref_allele";
                ($t, $p) = ($alt_allele =~ /(\S+)\](\S+)\]/);
            } elsif ($alt_allele =~ /^\]/) {
                $p_relative_strand = "+";
                $p_relative_position = "before_ref_allele";
                ($p, $t) = ($alt_allele =~ /\](\S+)\](\S+)/);
            } elsif ($alt_allele =~ /^\[/) {
                $p_relative_strand = "-";
                $p_relative_position = "before_ref_allele";
                ($p, $t) = ($alt_allele =~ /\[(\S+)\[(\S+)/);
            } else {
                print "Unexpected ALT field in the input vcf file for BND definiation:\n";
                print "$_\n";
                print "Exit!\n";
                die;
            }
	    
            my ($mate_breakend_chr, $mate_breakend_start) = ($p =~ /([^:]+):([^:]+)/);
	    my $ref_start_adj;
	    my $mate_breakend_start_adj;
	    if (($p_relative_strand eq "+") and ($p_relative_position eq "after_ref_allele")) {
		$ref_start_adj = $ref_start + 1;
		$mate_breakend_start_adj = $mate_breakend_start - 1;
	    } elsif (($p_relative_strand eq "+") and ($p_relative_position eq "before_ref_allele")) {
		$ref_start_adj = $ref_start - 1;
		$mate_breakend_start_adj = $mate_breakend_start + 1;
	    } elsif (($p_relative_strand eq "-") and ($p_relative_position eq "before_ref_allele")) {
		$ref_start_adj = $ref_start - 1;
		$mate_breakend_start_adj = $mate_breakend_start - 1;
	    } else {
		# (($p_relative_strand eq "-") and ($p_relative_position eq "after_ref_allele")) 
		$ref_start_adj = $ref_start + 1;
		$mate_breakend_start_adj = $mate_breakend_start + 1;
	    }

	    my @sv_breakend = sort ("$ref_chr:$ref_start", "$mate_breakend_chr:$mate_breakend_start", "$ref_chr:$ref_start_adj", "$mate_breakend_chr:$mate_breakend_start_adj");
	    foreach my $breakend (@sv_breakend) {
		my ($breakend_chr, $breakend_start) = split /:/, $breakend;
		if (not exists $sv{'breakend'}{$breakend_chr}{$breakend_start}{'sv_event_type_hash'}{$sv_event_type}) {
		    $sv{'breakend'}{$breakend_chr}{$breakend_start}{'sv_event_type_hash'}{$sv_event_type} = 1;
		}
		if (not exists $sv{'breakend'}{$breakend_chr}{$breakend_start}{'sv_event_name_hash'}{$sv_event_name}) {
		    $sv{'breakend'}{$breakend_chr}{$breakend_start}{'sv_event_name_hash'}{$sv_event_name} = 1;
		}
		foreach my $mate_breakend (@sv_breakend) {
		    if (not exists $sv{'breakend'}{$breakend_chr}{$breakend_start}{'sv_event_mate_breakend_hash'}{$mate_breakend}) {
			$sv{'breakend'}{$breakend_chr}{$breakend_start}{'sv_event_mate_breakend_hash'}{$mate_breakend} = 1;
		    }
		}
	    }
	}
    }

    # creating breakend group (nearby breakends within the merging distance)
    %{$sv{'breakend_grp'}} = ();
    if (keys %{$sv{'breakend'}}) {
	my $grp_index = 0;
	foreach my $breakend_chr (sort keys %{$sv{'breakend'}}) {
	    $grp_index++;
	    my @breakend_start = sort {$a <=> $b} keys %{$sv{'breakend'}{$breakend_chr}};
	    my $last_breakend_start = "";
	    foreach my $breakend_start (@breakend_start) {
		if ($last_breakend_start ne "") {
		    my $d = abs($breakend_start - $last_breakend_start);
		    if ($d > $max_offset) {
			$grp_index++;
		    }
		}
		my $breakend_grp_id = "breakend_group_${grp_index}";
		$sv{'breakend_grp'}{$breakend_grp_id}{'breakend_hash'}{$breakend_chr}{$breakend_start} = 1;
		$sv{'breakend'}{$breakend_chr}{$breakend_start}{'breakend_grp_id'} = $breakend_grp_id;
		$last_breakend_start = $breakend_start;
	    }
	}
    }


    # creating sv_event_mate_breakend_grp_mate_map
    my %sv_event_mate_breakend_grp_mate_map  = ();
    foreach my $breakend_grp_id (sort keys %{$sv{'breakend_grp'}}) {
	my %involved_breakend_grp_id = ();
	$involved_breakend_grp_id{$breakend_grp_id} = 1;
	foreach my $breakend_chr (sort keys %{$sv{'breakend_grp'}{$breakend_grp_id}{'breakend_hash'}}) {
	    foreach my $breakend_start (sort {$a <=> $b} keys %{$sv{'breakend_grp'}{$breakend_grp_id}{'breakend_hash'}{$breakend_chr}}) {
		my @sv_event_mate_breakend = sort keys %{$sv{'breakend'}{$breakend_chr}{$breakend_start}{'sv_event_mate_breakend_hash'}};
		foreach my $sv_event_mate_breakend (@sv_event_mate_breakend) {
		    my ($sv_event_mate_breakend_chr, $sv_event_mate_breakend_start) = split /:/, $sv_event_mate_breakend;
		    my $sv_event_mate_breakend_grp_id = $sv{'breakend'}{$sv_event_mate_breakend_chr}{$sv_event_mate_breakend_start}{'breakend_grp_id'};
		    if (not exists $involved_breakend_grp_id{$sv_event_mate_breakend_grp_id}) {
			$involved_breakend_grp_id{$sv_event_mate_breakend_grp_id} = 1;
		    }
		}
	    }
	}
	
	foreach my $previously_examined_breakend_grp_id (sort keys %sv_event_mate_breakend_grp_mate_map) {
	    if (exists $involved_breakend_grp_id{$previously_examined_breakend_grp_id}) {
		foreach my $previously_examined_mate_breakend_grp_id (sort keys %{$sv_event_mate_breakend_grp_mate_map{$previously_examined_breakend_grp_id}}) {
		    if (not exists $involved_breakend_grp_id{$previously_examined_mate_breakend_grp_id}) {
			$involved_breakend_grp_id{$previously_examined_mate_breakend_grp_id} = 1;
		    }
		}
	    }
	}
	
	foreach my $involved_breakend_grp_id (sort keys %involved_breakend_grp_id) {
	    foreach my $involved_breakend_grp_id2 (sort keys %involved_breakend_grp_id) {
		if (not exists $sv_event_mate_breakend_grp_mate_map{$involved_breakend_grp_id}{$involved_breakend_grp_id2}) {
		    $sv_event_mate_breakend_grp_mate_map{$involved_breakend_grp_id}{$involved_breakend_grp_id2} = 1
		}
	    }
	}
    }

    my %sv_event_type = ();
    my %sv_event_name = ();
    my %sv_event_mate_breakend = ();
    
    foreach my $breakend_grp_id (sort keys %{$sv{'breakend_grp'}}) {
	foreach my $sv_event_mate_breakend_grp_id (sort keys %{$sv_event_mate_breakend_grp_mate_map{$breakend_grp_id}}) {
	    foreach my $sv_event_mate_breakend_chr (sort keys %{$sv{'breakend_grp'}{$sv_event_mate_breakend_grp_id}{'breakend_hash'}}) {
		foreach my $sv_event_mate_breakend_start (sort {$a <=> $b} keys %{$sv{'breakend_grp'}{$sv_event_mate_breakend_grp_id}{'breakend_hash'}{$sv_event_mate_breakend_chr}}) {
		    my @sv_event_mate_breakend_sv_event_type = sort keys %{$sv{'breakend'}{$sv_event_mate_breakend_chr}{$sv_event_mate_breakend_start}{'sv_event_type_hash'}};
		    foreach my $sv_event_mate_breakend_sv_event_type (@sv_event_mate_breakend_sv_event_type) {
			$sv_event_type{$breakend_grp_id}{$sv_event_mate_breakend_sv_event_type} = 1;
		    }
		    my @sv_event_mate_breakend_sv_event_name = sort keys %{$sv{'breakend'}{$sv_event_mate_breakend_chr}{$sv_event_mate_breakend_start}{'sv_event_name_hash'}};
		    foreach my $sv_event_mate_breakend_sv_event_name (@sv_event_mate_breakend_sv_event_name) {
			$sv_event_name{$breakend_grp_id}{$sv_event_mate_breakend_sv_event_name} = 1;
		    }
		    if (not exists $sv_event_mate_breakend{$breakend_grp_id}{"$sv_event_mate_breakend_chr:$sv_event_mate_breakend_start"}) {
			$sv_event_mate_breakend{$breakend_grp_id}{"$sv_event_mate_breakend_chr:$sv_event_mate_breakend_start"} = 1;
		    }
		}
	    }
	}
    }
    foreach my $breakend_grp_id (sort keys %{$sv{'breakend_grp'}}) {
	foreach my $breakend_chr (sort keys %{$sv{'breakend_grp'}{$breakend_grp_id}{'breakend_hash'}}) {
	    foreach my $breakend_start (sort {$a <=> $b} keys %{$sv{'breakend_grp'}{$breakend_grp_id}{'breakend_hash'}{$breakend_chr}}) {
		%{$sv{'breakend'}{$breakend_chr}{$breakend_start}{'sv_event_mate_breakend_hash'}} = %{$sv_event_mate_breakend{$breakend_grp_id}}; 
		%{$sv{'breakend'}{$breakend_chr}{$breakend_start}{'sv_event_type_hash'}} = %{$sv_event_type{$breakend_grp_id}}; 
		%{$sv{'breakend'}{$breakend_chr}{$breakend_start}{'sv_event_name_hash'}} = %{$sv_event_name{$breakend_grp_id}};
	    } 
	}
	%{$sv{'breakend_grp'}{$breakend_grp_id}{'sv_event_mate_breakend_grp_hash'}} = %{$sv_event_mate_breakend_grp_mate_map{$breakend_grp_id}};
    }


    # create global sv_event_id
    if (keys %{$sv{'breakend'}}) {
	foreach my $breakend_chr (sort keys %{$sv{'breakend'}}) {
	    foreach my $breakend_start (sort {$a <=> $b} keys %{$sv{'breakend'}{$breakend_chr}}) {
		my $sv_event_mate_breakend = join ";", (sort keys %{$sv{'breakend'}{$breakend_chr}{$breakend_start}{'sv_event_mate_breakend_hash'}});
		my $sv_event_id  = $sv_event_mate_breakend;
		my $sv_event_type = join ";", (sort keys %{$sv{'breakend'}{$breakend_chr}{$breakend_start}{'sv_event_type_hash'}});
		my $sv_event_name = join ";", (sort keys %{$sv{'breakend'}{$breakend_chr}{$breakend_start}{'sv_event_name_hash'}});
		my $breakend_grp_id = $sv{'breakend'}{$breakend_chr}{$breakend_start}{'breakend_grp_id'};
		$sv{'breakend'}{$breakend_chr}{$breakend_start}{'sv_event_mate_breakend'} = $sv_event_mate_breakend;
		$sv{'breakend'}{$breakend_chr}{$breakend_start}{'sv_event_id'} = $sv_event_id;
		$sv{'breakend'}{$breakend_chr}{$breakend_start}{'sv_event_type'} = $sv_event_type;
		$sv{'breakend'}{$breakend_chr}{$breakend_start}{'sv_event_name'} = $sv_event_name;
		$sv{'breakend_grp'}{$breakend_grp_id}{'sv_event_mate_breakend'} = $sv_event_mate_breakend;
		$sv{'breakend_grp'}{$breakend_grp_id}{'sv_event_id'} = $sv_event_id;
		$sv{'breakend_grp'}{$breakend_grp_id}{'sv_event_type'} = $sv_event_type;
		$sv{'breakend_grp'}{$breakend_grp_id}{'sv_event_name'} = $sv_event_name;
	    }
	}
    }
    
    # print Dumper(%{$sv{'breakend'}}), "\n";

    %{$sv{'event'}} = ();
    if (keys %{$sv{'breakend_grp'}}) {
	foreach my $breakend_grp_id (sort keys %{$sv{'breakend_grp'}}) {
	    my $sv_event_id = $sv{'breakend_grp'}{$breakend_grp_id}{'sv_event_id'};
	    my $sv_event_type = $sv{'breakend_grp'}{$breakend_grp_id}{'sv_event_type'};
	    my $sv_event_name = $sv{'breakend_grp'}{$breakend_grp_id}{'sv_event_name'};
	    if (not exists $sv{'event'}{$sv_event_id}) {
                $sv{'event'}{$sv_event_id}{'sv_event_id'} = $sv_event_id;
                $sv{'event'}{$sv_event_id}{'sv_event_name'} = $sv_event_name;
                $sv{'event'}{$sv_event_id}{'sv_event_type'} = $sv_event_type;
	    }
	    $sv{'event'}{$sv_event_id}{'sv_event_mate_breakend_grp_hash'}{$breakend_grp_id} = 1;
	}
	foreach my $sv_event_id (sort keys %{$sv{'event'}}) {
	    foreach my $sv_event_breakend_grp_id (sort keys %{$sv{'event'}{$sv_event_id}{'sv_event_mate_breakend_grp_hash'}}) {
		foreach my $sv_event_breakend_grp_chr (sort keys %{$sv{'breakend_grp'}{$sv_event_breakend_grp_id}{'breakend_hash'}}) {
		    my @sv_event_breakend_grp_start = sort {$a <=> $b} keys %{$sv{'breakend_grp'}{$sv_event_breakend_grp_id}{'breakend_hash'}{$sv_event_breakend_grp_chr}};
		    my $sv_event_breakend_grp_mean_start = cal_mean(\@sv_event_breakend_grp_start);
		    $sv{'event'}{$sv_event_id}{'sv_event_mate_breakend_grp_rep_hash'}{$sv_event_breakend_grp_id}{"$sv_event_breakend_grp_chr:$sv_event_breakend_grp_mean_start"} = 1;
		}
	    }
	}
    }
    return %sv;
}

sub cal_mean {
    my $data_arrayref = shift @_;
    my $data_num = scalar @$data_arrayref;
    my $data_sum = sum @$data_arrayref;
    my $mean;
    if ($data_num == 0) {
	$mean = "NA";
    } else {
	$mean = $data_sum/$data_num;
    }
    return $mean;
}

sub count_variants_for_SV_vcf {
    my $vcf_hashref = shift @_;
    my %n = ();
    $n{'event'} = scalar (keys %{$$vcf_hashref{'event'}});
    $n{'breakend_grp'} = scalar (keys %{$$vcf_hashref{'breakend_grp'}});
    $n{'breakend'} = 0;
    foreach my $chr (sort keys %{$$vcf_hashref{'breakend'}}) {
	my $count = scalar (keys %{$$vcf_hashref{'breakend'}{$chr}});
	$n{'breakend'} += $count;
    }
    return %n;
}

sub SV_vcf_intersect {
    my ($ref_vcf_hashref, $query_vcf_hashref, $max_offset) = @_;
    my %intersect_vcf = ();
    foreach my $ref_chr (sort keys %{$$ref_vcf_hashref{'breakend'}}) {
	foreach my $ref_start (sort {$a <=> $b} keys %{$$ref_vcf_hashref{'breakend'}{$ref_chr}}) {
	    if (exists $$query_vcf_hashref{'breakend'}{$ref_chr}) {
		foreach my $query_start (sort {$a <=> $b} keys %{$$query_vcf_hashref{'breakend'}{$ref_chr}}) {
		    my $d = abs($query_start - $ref_start);
		    if ($d <= $max_offset) {
			$intersect_vcf{'breakend'}{$ref_chr}{$ref_start}{'match_query_breakend_hash'}{$ref_chr}{$query_start} = $d;
			my $match_query_breakend_grp_id = $$query_vcf_hashref{'breakend'}{$ref_chr}{$query_start}{'breakend_grp_id'};
			$intersect_vcf{'breakend'}{$ref_chr}{$ref_start}{'match_query_breakend_grp_hash'}{$match_query_breakend_grp_id}{'breakend_hash'}{$ref_chr}{$query_start} = $d;
		    }
		}
	    }
	}
    }
    
    my %match_query_breakend_grp_id = ();
    foreach my $ref_breakend_grp_id (sort keys %{$$ref_vcf_hashref{'breakend_grp'}}) {
	foreach my $ref_chr (sort keys %{$$ref_vcf_hashref{'breakend_grp'}{$ref_breakend_grp_id}{'breakend_hash'}}) {
	    foreach my $ref_start (sort {$a <=> $b} keys %{$$ref_vcf_hashref{'breakend_grp'}{$ref_breakend_grp_id}{'breakend_hash'}{$ref_chr}}) {
		if (exists $intersect_vcf{'breakend'}{$ref_chr}{$ref_start}) {
		    my @match_query_breakend_grp_id = sort keys %{$intersect_vcf{'breakend'}{$ref_chr}{$ref_start}{'match_query_breakend_grp_hash'}};
		    foreach my $match_query_breakend_grp_id (@match_query_breakend_grp_id) {
			if (not exists $intersect_vcf{'breakend_grp'}{$ref_breakend_grp_id}{$match_query_breakend_grp_id}) {
			    $intersect_vcf{'breakend_grp'}{$ref_breakend_grp_id}{$match_query_breakend_grp_id} = 1;
			    $match_query_breakend_grp_id{$match_query_breakend_grp_id}{$ref_breakend_grp_id} = 1;
			}
		    }
		}
	    }
	}
    }
    
    foreach my $ref_event_id (sort keys %{$$ref_vcf_hashref{'event'}}) {
	my $ref_based_check_flag = 0;
	my $query_based_check_flag = 0;
	my %match_query_event = ();
	my @ref_event_breakend_grp_id = sort keys %{$$ref_vcf_hashref{'event'}{$ref_event_id}{'sv_event_mate_breakend_grp_hash'}};
	my %ref_event_breakend_grp_id = ();
	foreach my $ref_event_breakend_grp_id (@ref_event_breakend_grp_id) {
	    $ref_event_breakend_grp_id{$ref_event_breakend_grp_id} = 1;
	    if (not exists $intersect_vcf{'breakend_grp'}{$ref_event_breakend_grp_id}) {
		$ref_based_check_flag = 1;
		last;
	    } else {
		my @match_query_breakend_grp_id = sort keys %{$intersect_vcf{'breakend_grp'}{$ref_event_breakend_grp_id}};
		foreach my $match_query_breakend_grp_id (@match_query_breakend_grp_id) {
		    foreach my $match_query_breakend_chr (sort keys %{$$query_vcf_hashref{'breakend_grp'}{$match_query_breakend_grp_id}{'breakend_hash'}}) {
			foreach my $match_query_breakend_start (sort {$a <=> $b} keys %{$$query_vcf_hashref{'breakend_grp'}{$match_query_breakend_grp_id}{'breakend_hash'}{$match_query_breakend_chr}}) {
			    my $match_query_event_id = $$query_vcf_hashref{'breakend'}{$match_query_breakend_chr}{$match_query_breakend_start}{'sv_event_id'};
			    if (not exists $match_query_event{$match_query_event_id}) {
				$match_query_event{$match_query_event_id}{$ref_event_id} = 1;
			    }
			}
		    }
		}
	    }
	}
	my @match_query_event_id = sort keys %match_query_event;
	my $match_query_event_id;
	if ((scalar @match_query_event_id) == 1) {
	    $match_query_event_id = shift @match_query_event_id;
	    foreach my $query_event_breakend_grp_id (sort keys %{$$query_vcf_hashref{'event'}{$match_query_event_id}{'sv_event_mate_breakend_grp_hash'}}) {
		if (not exists $match_query_breakend_grp_id{$query_event_breakend_grp_id}) {
		    $query_based_check_flag = 1;
		    last;
		} else {
		    my @match_ref_breakend_grp_id = sort keys %{$match_query_breakend_grp_id{$query_event_breakend_grp_id}};
		    foreach my $match_ref_breakend_grp_id (@match_ref_breakend_grp_id) {
			if (not exists $ref_event_breakend_grp_id{$match_ref_breakend_grp_id}) {
			    $query_based_check_flag = 1;
			    last;
			}
		    }
		}
	    }
	} else {
	    $query_based_check_flag = 1;
	}
	
	if (($ref_based_check_flag == 0) and ($query_based_check_flag == 0)) {
	    $intersect_vcf{'event'}{$ref_event_id}{'sv_event_mate_breakend_hash'} = $$ref_vcf_hashref{'event'}{$ref_event_id}{'sv_event_mate_breakend_hash'};
	    $intersect_vcf{'event'}{$ref_event_id}{'sv_event_mate_breakend_grp_hash'} = $$ref_vcf_hashref{'event'}{$ref_event_id}{'sv_event_mate_breakend_grp_hash'};
	    $intersect_vcf{'event'}{$ref_event_id}{'match_query_event_id'} = $match_query_event_id;
	}
    }
    return %intersect_vcf;
}


sub SV_vcf_uniq {
    my ($ref_vcf_hashref, $query_vcf_hashref, $max_offset) = @_;
    my %uniq_vcf = ();
    foreach my $ref_chr (sort keys %{$$ref_vcf_hashref{'breakend'}}) {
	foreach my $ref_start (sort {$a <=> $b} keys %{$$ref_vcf_hashref{'breakend'}{$ref_chr}}) {
	    $uniq_vcf{'breakend'}{$ref_chr}{$ref_start} = 1;
	    $uniq_vcf{'breakend'}{$ref_chr}{$ref_start} = 1;
	    if (exists $$query_vcf_hashref{'breakend'}{$ref_chr}) {
		foreach my $query_start (sort {$a <=> $b} keys %{$$query_vcf_hashref{'breakend'}{$ref_chr}}) {
		    my $d_tmp = abs($query_start - $ref_start);
		    if ($d_tmp <= $max_offset) {
			delete $uniq_vcf{'breakend'}{$ref_chr}{$ref_start};
		    } 
		}
	    }
	}
    }
    
    foreach my $ref_breakend_grp_id (sort keys %{$$ref_vcf_hashref{'breakend_grp'}}) {
	my $flag = 0;
	foreach my $ref_chr (sort keys %{$$ref_vcf_hashref{'breakend_grp'}{$ref_breakend_grp_id}{'breakend_hash'}}) {
	    foreach my $ref_start (sort {$a <=> $b} keys %{$$ref_vcf_hashref{'breakend_grp'}{$ref_breakend_grp_id}{'breakend_hash'}{$ref_chr}}) {
		if (not exists $uniq_vcf{'breakend'}{$ref_chr}{$ref_start}) {
		    $flag = 1;
		    last;
		}
	    }
	}
	if ($flag == 0) {
	    foreach my $ref_chr (sort keys %{$$ref_vcf_hashref{'breakend_grp'}{$ref_breakend_grp_id}{'breakend_hash'}}) {
		foreach my $ref_start (sort {$a <=> $b} keys %{$$ref_vcf_hashref{'breakend_grp'}{$ref_breakend_grp_id}{'breakend_hash'}{$ref_chr}}) {
		    $uniq_vcf{'breakend_grp'}{$ref_breakend_grp_id}{'breakend_hash'}{$ref_chr}{$ref_start} = 1;
		}
	    }
	}
    }
    return %uniq_vcf;
}



#-----------------------------------------------------------------
#----------------  Documentation / Usage / Help ------------------

=head1 NAME

vcf_benchmarker.pl - compare the reference and query vcf files to determine the precision, recall, and F1 score of the query variant calling result by assuming the reference vcf file as the ground truth.

=head1 SYNOPSIS

perl vcf_benchmarker.pl [options] [file ...]

=head1 OPTIONS

=over 8

=item B<-help> or B<-h>

Print help message. Example: -h.

=item B<-man> or B<-m>

Print more detailed help message. Example: -m.

=item B<-version> or B<-v>

Print version information. Example: -v.

=item B<-ref_vcf> or B<-r> 

Specify the reference variant calling file (in vcf or vcf.gz format) to be used as the ground truth. Example: -ref_vcf ref.vcf(.gz).

=item B<-query_vcf> or B<-q>

Specify the query variant calling file (in vcf or vcf.gz format) to be used for comparing with the the reference vcf. Example: -query_vcf query.vcf(.gz).

=item B<-vcf_for_sv>

Specify whether the input vcf files are used for denoting structural variants (SVs). Default: '-vcf_for_sv no'.

=item B<-snp_indel_qual_filter>

The quality score filter for SNP/INDEL records in vcf files. Default: '-snp_indel_qual_filter 30'.

=item B<-sv_qual_filter>

The quality filter for SV records in vcf files. Mulitple filters separated by ";" can be defined in the same time. Default: -sv_qual_filter LowQual;q5'.

=item B<-max_offset> or B<-d>

The maxium allowance for the breakend coordinate offset in basepair when comparing reference and query SV vcf files. SVs with breakends off by more than this specified distance will be considered different. Default: '-max_offset 10' (i.e. 10 bp).

=item B<-prefix> or B<-p>

Specify the file name prefix for the output files. Example: -prefix test_prefix. Default: '-prefix output_prefix'.

=back

=head1 DESCRIPTION

B<vcf_benchmarker.pl> compare the reference and query vcf files to determine the precision, recall, and F1 score of the query variant calling result by assuming the reference vcf file as the ground truth.
    
=head1 AUTHOR

B<Jia-Xing Yue> (GitHub ID: yjx1217)                                                              

=head1 VERSION

B<version> v2019.05.09

=cut
