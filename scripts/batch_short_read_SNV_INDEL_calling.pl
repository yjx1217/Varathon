#!/usr/bin/perl
use warnings;
use strict;
use Getopt::Long;
use Env;
use Cwd;

##############################################################
#  script: batch_short_read_SNV_INDEL_calling.pl
#  author: Jia-Xing Yue (GitHub ID: yjx1217)
#  last edited: 2023.05.08
#  description: run short-read-based SNV & INDEL calling for a batch of samples
#  example: perl batch_short_read_SNV_INDEL_calling.pl -i Master_Sample_Table.txt -threads 4 -b $batch_id -ref_genome ref_genome.fa  -short_read_mapping_dir ./../01.Short_Read_Mapping -min_mapping_quality 30 -min_variant_calling_quality 30 -caller GATK4  -ploidy 1 -excluded_chr_list yeast.excluded_chr_list.txt
##############################################################

my $VARATHON_HOME = $ENV{VARATHON_HOME};
my $samtools_dir = $ENV{samtools_dir};
my $picard_dir = $ENV{picard_dir};
my $clair3_dir = $ENV{clair3_dir};
my $gatk4_dir = $ENV{gatk4_dir};
my $freebayes_dir = $ENV{freebayes_dir};
my $parallel_dir = $ENV{parallel_dir};
my $vt_dir = $ENV{vt_dir};
my $vcflib_dir = $ENV{vcflib_dir};
my $tabix_dir = $ENV{tabix_dir};
my $sample_table = "Master_Sample_Table.txt";
my $batch_id;
my $threads = 1;
my $ref_genome;
my $short_read_mapping_dir;
my $min_mapping_quality = 30;
my $min_variant_calling_quality = 30;
my $ploidy = 2;
my $caller = "gatk4"; # "freebayes"
my $debug = "no";
my $excluded_chr_list;

GetOptions('sample_table|i:s' => \$sample_table,
	   'threads|t:i' => \$threads,
	   'batch_id|b:s' => \$batch_id,
	   'ref_genome|ref_genome:s' => \$ref_genome,
	   'short_read_mapping_dir|short_read_mapping_dir:s' => \$short_read_mapping_dir,
	   'min_mapping_quality|mq:i' => \$min_mapping_quality,
	   'min_variant_calling_quality|vq:f' => \$min_variant_calling_quality,
	   'caller|c:s' => \$caller,
	   'ploidy|p:i' => \$ploidy,
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
    my $local_time = localtime();
    print "[$local_time] processing sample $sample_id with read-alignment filtering\n";
    my $ref_genome_file = "$base_dir/$ref_genome";
    my $bam_file = "$base_dir/$short_read_mapping_dir/$batch_id/$sample_id/$sample_id.final.bam";
    print "Check the specified reference genome file:\n";
    if (-e $ref_genome_file) {
        print "Successfully located the specified reference genome file: $ref_genome_file.\n";
    } else {
        print "Cannot find the specified reference genome file: $ref_genome_file!\n";
        print "Exit!\n";
        exit;
    }
    print "Check the specified short-read mapping bam file:\n";
    if (-e $bam_file) {
        print "Successfully located the specified short-read mapping bam file: $bam_file.\n";
    } else {
        print "Cannot find the specified short-read mapping bam file: $bam_file!\n";
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
    ## filter bam file by mapping quality
    system("$samtools_dir/samtools view -b -q $min_mapping_quality $base_dir/$short_read_mapping_dir/$batch_id/$sample_id/$sample_id.final.bam  >$sample_id.filtered.bam");
    # index the filtered.bam file
    system("$samtools_dir/samtools index $sample_id.filtered.bam");

    $local_time = localtime();
    print "[$local_time] processing sample $sample_id with SNV and INDEL calling\n";

    # SNV and INDEL calling 
    if ($caller eq "gatk4") {
	system("/usr/bin/time -v java -Djava.io.tmpdir=./tmp -XX:ParallelGCThreads=$threads -jar $gatk4_dir/gatk4.jar HaplotypeCaller -R ref.genome.fa -I $sample_id.filtered.bam -O $sample_id.$caller.raw.vcf -ploidy $ploidy");
	system("rm $sample_id.$caller.raw.vcf.idx");
    } elsif ($caller eq "freebayes") {
	system("/usr/bin/time -v $freebayes_dir/python $freebayes_dir/fasta_generate_regions.py ref.genome.fa.fai 100000 > $sample_id.region.txt");
	system("cat $sample_id.region.txt | /usr/bin/time -v $parallel_dir/parallel -k -j $threads $freebayes_dir/freebayes -f ref.genome.fa  -p $ploidy $sample_id.filtered.bam --region {} | $freebayes_dir/python $freebayes_dir/vcffirstheader | $freebayes_dir/vcfstreamsort -w 1000 | $freebayes_dir/vcfuniq > $sample_id.$caller.raw.vcf");
    } elsif ($caller eq "clair3") {
	if ($ploidy eq "2") {
	    system("/usr/bin/time -v $clair3_dir/run_clair3.sh --threads=${threads} --model_path=$clair3_dir/models/ilmn --platform=ilmn --ref_fn=ref.genome.fa --include_all_ctgs --bam_fn=$sample_id.filtered.bam --sample_name=$sample_id --output=$sample_id.$caller.output_dir");
	} else {
	    system("/usr/bin/time -v $clair3_dir/run_clair3.sh --threads=${threads} --model_path=$clair3_dir/models/ilmn --platform=ilmn --ref_fn=ref.genome.fa --include_all_ctgs --no_phasing_for_fa --bam_fn=$sample_id.filtered.bam --sample_name=$sample_id --output=$sample_id.$caller.output_dir");
	}
	system("/usr/bin/time -v gzip -dc $sample_id.$caller.output_dir/merge_output.vcf.gz | $vcflib_dir/vcfstreamsort | $vt_dir/vt sort - > $sample_id.$caller.raw.vcf");
	if ($debug eq "no") {
	    system("rm -r $sample_id.$caller.output_dir");
	}
    } else {
	die "Unknown caller: $caller! Please only use \"gatk4\" or \"freebayes\" or \"clair3\" as your short-read-based SNV & INDEL caller!\n";
    }

    $local_time = localtime();
    print "[$local_time] processing sample $sample_id with variant calls filtering\n";

    system("$vt_dir/vt decompose_blocksub $sample_id.$caller.raw.vcf -a -o $sample_id.$caller.decomposed.vcf");
    system("$vt_dir/vt normalize $sample_id.$caller.decomposed.vcf -r ref.genome.fa -f \"~VARIANT_CONTAINS_N\"| $vt_dir/vt uniq - -o  $sample_id.$caller.normalized.vcf");

    system("$vt_dir/vt annotate_variants $sample_id.$caller.normalized.vcf -r ref.genome.fa -o $sample_id.$caller.vt_annotated.vcf");
    system("$vt_dir/vt view $sample_id.$caller.vt_annotated.vcf -f \"VTYPE==SNP\" -o $sample_id.$caller.vt_annotated.SNV.vcf");
    system("$vt_dir/vt view $sample_id.$caller.vt_annotated.vcf -f \"VTYPE==INDEL\" -o $sample_id.$caller.vt_annotated.INDEL.vcf");
    system("$vcflib_dir/vcffilter -f \"QUAL > $min_variant_calling_quality\" $sample_id.$caller.vt_annotated.vcf > $sample_id.$caller.filtered.vcf");
    system("$vcflib_dir/vcffilter -f \"QUAL > $min_variant_calling_quality\" $sample_id.$caller.vt_annotated.SNV.vcf > $sample_id.$caller.filtered.SNV.vcf");
    system("$vcflib_dir/vcffilter -f \"QUAL > $min_variant_calling_quality\" $sample_id.$caller.vt_annotated.INDEL.vcf > $sample_id.$caller.filtered.INDEL.vcf");

    # system("$tabix_dir/bgzip -c $sample_id.$caller.filtered.SNV.vcf > $sample_id.$caller.filtered.SNV.vcf.gz");
    # system("$tabix_dir/tabix -p vcf $sample_id.$caller.filtered.SNV.vcf.gz");
    # system("$tabix_dir/bgzip -c $sample_id.$caller.filtered.INDEL.vcf > $sample_id.$caller.filtered.INDEL.vcf.gz");
    # system("$tabix_dir/tabix -p vcf $sample_id.$caller.filtered.INDEL.vcf.gz");

    $local_time = localtime();
    print "[$local_time] finishing short-read-based SNV and INDEL calling for sample $sample_id \n";
  
    # remove intermediate files
    if ($debug eq "no") {
	system("rm ref.genome.fa");
	system("rm ref.genome.fa.fai");
	system("rm ref.genome.dict");
	system("rm $sample_id.filtered.bam");
	system("rm $sample_id.filtered.bam.bai");
	system("rm $sample_id.$caller.raw.vcf");
	system("rm $sample_id.$caller.decomposed.vcf");
	system("rm $sample_id.$caller.normalized.vcf");
	# system("rm $sample_id.$caller.vt_annotated.vcf");
	# system("rm $sample_id.$caller.vt_annotated.SNV.vcf");
	# system("rm $sample_id.$caller.vt_annotated.INDEL.vcf");
	# system("rm $sample_id.$caller.filtered.SNV.vcf");
	# system("rm $sample_id.$caller.filtered.INDEL.vcf");
	system("rm -r tmp");
	if ($caller eq "freebayes") {
	    system("rm $sample_id.region.txt")
	}
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
        my ($sample_id, $PE_read_files, $note) = split /\s+/, $_;
	push @$sample_table_arrayref, $sample_id;
        my ($R1_read_file, $R2_read_file) = split /,/, $PE_read_files;
        $$sample_table_hashref{$sample_id}{'R1_read_file'} = $R1_read_file;
        $$sample_table_hashref{$sample_id}{'R2_read_file'} = $R2_read_file;
	$$sample_table_hashref{$sample_id}{'note'} = $note;
    }
}

