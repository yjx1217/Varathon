#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Env;
use Cwd;

##############################################################
#  script: batch_long_read_SV_calling.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2019.12.28
#  description: run long-read-based SV calling for a batch of samples
#  example: perl batch_long_read_SV_calling.pl -i Master_Sample_Table.txt -threads 4 -b $batch_id -ref_genome ref_genome.fa  -long_read_mapping_dir ./../01.Long_Read_Mapping -min_mapping_quality 30 -caller sniffles -excluded_chr_list yeast.excluded_chr_list.txt -long_read_technology pacbio
##############################################################

my $VARATHON_HOME = $ENV{VARATHON_HOME};
my $samtools_dir = $ENV{samtools_dir};
my $picard_dir = $ENV{picard_dir};
my $windowmasker_dir = $ENV{windowmasker_dir};
my $bedtools_dir = $ENV{bedtools_dir};
my $sniffles_dir = $ENV{sniffles_dir};
my $svim_dir = $ENV{svim_dir};
my $picky_dir = $ENV{picky_dir};
my $pbsv_dir = $ENV{pbsv_dir};
my $nanosv_dir = $ENV{nanosv_dir};
my $vt_dir = $ENV{vt_dir};
my $vcflib_dir = $ENV{vcflib_dir};
my $tabix_dir = $ENV{tabix_dir};
my $sample_table = "Master_Sample_Table.txt";
my $batch_id;
my $threads = 1;
my $ref_genome;
my $long_read_technology = "pacbio"; 
my $long_read_mapping_dir;
my $min_mapping_quality = 30;
my $min_read_support = 5;
# my $ploidy = 2;
my $caller = "sniffles"; # "sniffles"
my $debug = "no";
my $excluded_chr_list;

GetOptions('sample_table|i:s' => \$sample_table,
	   'threads|t:i' => \$threads,
	   'batch_id|b:s' => \$batch_id,
	   'ref_genome|ref_genome:s' => \$ref_genome,
	   'long_read_technology|long_read_technology:s' => \$long_read_technology,
	   'long_read_mapping_dir|long_read_mapping_dir:s' => \$long_read_mapping_dir,
	   'min_mapping_quality|mq:i' => \$min_mapping_quality,
	   'min_read_support|mr:i' => \$min_read_support,
	   'caller|c:s' => \$caller,
	   'debug|d:s' => \$debug,
	   'excluded_chr_list|excluded_chr_list:s' => \$excluded_chr_list); 

my $sample_table_fh = read_file($sample_table);
my %sample_table = ();
my @sample_table = ();
parse_sample_table($sample_table_fh, \%sample_table, \@sample_table);
my $base_dir = cwd();
my $output_dir = "$batch_id";
system("mkdir $output_dir");

my $sample_count = 0;
foreach my $sample_id (@sample_table) {
    $sample_count++;
    print "\n";
    my $local_time = localtime();
    print "[$local_time] processing sample $sample_id with read-alignment filtering\n";
    my $ref_genome_file = "$base_dir/$ref_genome";
    my $bam_file = "$base_dir/$long_read_mapping_dir/$batch_id/$sample_id/$sample_id.dedup.bam";
    print "Check the specified reference genome file:\n";
    if (-e $ref_genome_file) {
        print "Successfully located the specified reference genome file: $ref_genome_file.\n";
    } else {
        print "Cannot find the specified reference genome file: $ref_genome_file!\n";
        print "Exit!\n";
        exit;
    }
    print "Check the specified long-read mapping bam file:\n";
    if (-e $bam_file) {
        print "Successfully located the specified long-read mapping bam file: $bam_file.\n";
    } else {
        print "Cannot find the specified long-read mapping bam file: $bam_file!\n";
        print "Exit!\n";
        exit;
    }
    my $sample_output_dir = "$base_dir/$output_dir/$sample_id";
    system("mkdir -p  $sample_output_dir");
    chdir("$sample_output_dir") or die "cannot change directory to: $!\n";
    system("mkdir tmp");
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
    system("java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar CreateSequenceDictionary -REFERENCE ref.genome.fa -OUTPUT ref.genome.dict");
    if ($caller =~ /(sniffles|svim|pbsv|nanosv)/) {
	## filter bam file by mapping quality
	if ($caller =~ /pbsv/) {
	    system("$samtools_dir/samtools view -h -q $min_mapping_quality $base_dir/$long_read_mapping_dir/$batch_id/$sample_id/$sample_id.dedup.bam | awk '\$6 !~ /H/{print}'| $samtools_dir/samtools view -h -b - >$sample_id.filtered.bam");
	} else {
	    system("$samtools_dir/samtools view -h -b -q $min_mapping_quality $base_dir/$long_read_mapping_dir/$batch_id/$sample_id/$sample_id.dedup.bam  >$sample_id.filtered.bam");
	}
	## add read group
	system("java -Djava.io.tmpdir=./tmp -Dpicard.useLegacyParser=false -XX:ParallelGCThreads=$threads -jar $picard_dir/picard.jar AddOrReplaceReadGroups -I $sample_id.filtered.bam -O $sample_id.rdgrp.bam -SORT_ORDER coordinate -RGID $sample_id -RGLB $sample_id -RGPL $long_read_technology -RGPU $sample_id -RGSM $sample_id -RGCN 'RGCN'");
	## calculate MD tag
	system("$samtools_dir/samtools calmd -b --threads $threads $sample_id.rdgrp.bam ref.genome.fa > $sample_id.calmd.bam");
	# index the input bam file
	system("$samtools_dir/samtools index $sample_id.calmd.bam");
    }
    $local_time = localtime();
    print "[$local_time] processing sample $sample_id with SV calling\n";
    # SV calling 
    if ($caller eq "sniffles") {
	system("/usr/bin/time -v $sniffles_dir/sniffles -t $threads -m $sample_id.calmd.bam --min_length 50 --min_support $min_read_support -v $sample_id.$caller.SV.vcf");
    } elsif ($caller eq "svim") {
	system("/usr/bin/time -v $svim_dir/svim alignment --min_sv_size 50 --max_sv_size 1000000 ${sample_id}_svim_out $sample_id.calmd.bam");
	system("cp ./${sample_id}_svim_out/final_results.vcf $sample_id.$caller.SV.vcf");
    } elsif ($caller eq "picky") {
	if (-e "$base_dir/$long_read_mapping_dir/$batch_id/$sample_id/$sample_id.maf.gz") {
	    system("/usr/bin/time -v gzip -dc $base_dir/$long_read_mapping_dir/$batch_id/$sample_id/$sample_id.maf.gz | $picky_dir/picky.pl selectRep --thread $threads --preload 10 1 > $sample_id.align");
	    system("/usr/bin/time -v $samtools_dir/samtools bam2fq -n -i $base_dir/$long_read_mapping_dir/$batch_id/$sample_id/$sample_id.dedup.bam > $sample_id.fastq");
	} elsif (-e "$base_dir/$long_read_mapping_dir/$batch_id/$sample_id/$sample_id.dedup.bam") {
	    system("/usr/bin/time -v $samtools_dir/samtools view -h -q $min_mapping_quality $base_dir/$long_read_mapping_dir/$batch_id/$sample_id/$sample_id.dedup.bam |$samtools_dir/samtools sort -n -@ $threads - -o $sample_id.sam");
	    system("/usr/bin/time -v $samtools_dir/samtools bam2fq -i $base_dir/$long_read_mapping_dir/$batch_id/$sample_id/$sample_id.dedup.bam > $sample_id.fastq");
	    system("/usr/bin/time -v $picky_dir/picky.pl sam2align < $sample_id.sam > $sample_id.align");
	    if ($debug eq "no") {
		system("rm $sample_id.sam");
	    }
	} else {
	    print "Error! Cannot find expected input file: $base_dir/$long_read_mapping_dir/$batch_id/$sample_id/$sample_id.maf or $base_dir/$long_read_mapping_dir/$batch_id/$sample_id/$sample_id.dedup.bam !\n";
	    print "Exit!\n";
	    die;
	}
	system("/usr/bin/time -v $picky_dir/picky.pl callSV --oprefix $sample_id --genome ref.genome.fa --fastq $sample_id.fastq  < $sample_id.align");
	system("/usr/bin/time -v $picky_dir/picky.pl xls2vcf --xls $sample_id.profile.DEL.xls --xls $sample_id.profile.INS.xls --xls $sample_id.profile.INDEL.xls --xls $sample_id.profile.INV.xls --xls $sample_id.profile.TTLC.xls --xls $sample_id.profile.TDSR.xls --xls $sample_id.profile.TDC.xls  --re $min_read_support > $sample_id.$caller.SV.vcf");
	if ($debug eq "no") {
	    system("rm $sample_id.align");
	    system("rm $sample_id.fastq");
	}
    } elsif ($caller eq "pbsv") {
	system("/usr/bin/time -v $pbsv_dir/pbsv discover $sample_id.calmd.bam $sample_id.svsig.gz");
	system("/usr/bin/time -v $pbsv_dir/pbsv call --num-threads  $threads --min-sv-length 50 --call-min-reads-all-samples $min_read_support --call-min-reads-one-sample $min_read_support ref.genome.fa $sample_id.svsig.gz $sample_id.$caller.SV.vcf");
	if ($debug eq "no") {
	    system("rm $sample_id.svsig.gz");
	}
    } elsif ($caller eq "nanosv") {
	system("/usr/bin/time -v $windowmasker_dir/windowmasker -mk_counts -in ref.genome.fa -out ref.genome.masking_library.ustat");
	system("/usr/bin/time -v $windowmasker_dir/windowmasker -ustat ref.genome.masking_library.ustat -in ref.genome.fa -out ref.genome.softmask.fa -outfmt fasta -dust true");
	system("perl $VARATHON_HOME/scripts/softmask2hardmask.pl -i ref.genome.softmask.fa -o ref.genome.hardmask.fa");
	system("perl $VARATHON_HOME/scripts/find_motif_in_genome.pl -i ref.genome.hardmask.fa -m $VARATHON_HOME/data/hardmask.motif.txt -p ref.genome.hardmask");
	my $masking_details_txt_fh = read_file("ref.genome.hardmask.masking_details.txt");
	my $masking_details_bed_fh = write_file("ref.genome.hardmask.masking_details.bed");
	while (<$masking_details_txt_fh>) {
	    chomp;
	    /^motif\tmatch_id\tref/ and next;
	    my ($motif, $match_id, $match_chr, $match_start, $match_end, $match_case) = split /\t/, $_;
	    $match_start = $match_start - 1;
	    print $masking_details_bed_fh "$match_chr\t$match_start\t$match_end\n";
	}
	# Choose 1M random spots in the genome
	system("$bedtools_dir/bedtools random -l 1 -seed 20190518 -n 1000000 -g ref.genome.fa.fai > ref.genome.random_sample.raw.bed");
	# Shuffle positions while excluding gaps/simple repeats
	system("$bedtools_dir/bedtools shuffle -excl ref.genome.hardmask.masking_details.bed -noOverlapping -i ref.genome.random_sample.raw.bed -g ref.genome.fa.fai > ref.genome.random_sample.shuffled.bed");
	system("/usr/bin/time -v $nanosv_dir/NanoSV -t $threads -b ref.genome.random_sample.shuffled.bed -s $samtools_dir/samtools -o $sample_id.$caller.SV.vcf $sample_id.calmd.bam");
    } else {
	print "Error!!! Unknown caller: $caller! Please use one of the following as your long-read-based SV caller!\n";
	print "\"sniffles\" (for both PacBio and Nanopore)\n";
	print "\"svim\" (for both PacBio and Nanopore)\n";
	print "\"picky\" (for both PacBio and Nanopore)\n";
	print "\"pbsv\" (for PacBio only)\n";
	print "\"nanosv\" (for Nanopore only)\n";
	print "Exit!\n";
	die;
    }
    $local_time = localtime();
    print "[$local_time] finishing long-read-based SV calling for sample $sample_id \n";

    # remove intermediate files
    if ($debug eq "no") {
	if (-e "ref.genome.fa") {
	    system("rm ref.genome.fa");
	}
	if (-e "ref.genome.fa.fai") {
	    system("rm ref.genome.fa.fai");
	}
	if (-e "ref.genome.dict") {
	    system("rm ref.genome.dict");
	}
	if ($caller =~ /(sniffles|svim|pbsv|nanosv)/) {
	    system("rm $sample_id.filtered.bam");
	    system("rm $sample_id.rdgrp.bam");
	    system("rm $sample_id.calmd.bam");
	    system("rm $sample_id.calmd.bam.bai");
	}
	if ($caller eq "nanosv") {
	    system("rm ref.genome.masking_library.ustat");
	    system("rm ref.genome.softmask.fa");
	    system("rm ref.genome.hardmask.fa");
	    system("rm ref.genome.hardmask.masking_summary.txt");
	    system("rm ref.genome.hardmask.masking_details.txt");
	    system("rm ref.genome.hardmask.masking_details.bed");
	    system("rm ref.genome.random_sample.raw.bed");
	    system("rm ref.genome.random_sample.shuffled.bed");
	}
	system("rm -r tmp");
    }
    chdir("./../") or die "cannot change directory to: $!\n";
}


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
        my ($sample_id, $long_read_file, $note) = split /\s+/, $_;
	push @$sample_table_arrayref, $sample_id;
        $$sample_table_hashref{$sample_id}{'long_read_file'} = $long_read_file;
	$$sample_table_hashref{$sample_id}{'note'} = $note;
    }
}


