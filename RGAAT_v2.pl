#!~/miniconda2/bin/perl
#Author:Liu Wanfei <liuwf@big.ac.cn>
#Description: This program can assemble and/or annotate genome for new genome and known genome upgrade using sequence alignment file (SAM or BAM format), sequence variant file (VCF format or five coloum table (tab-delimited, including chromosome, position, id, reference allele and alternative allele)) or new genome sequence file (FASTA format) based on reference genome sequence file (FASTA format) and/or annotation file (TBL, GTF, GFF, GFF3 or BED format).

use strict;
#use warnings;
use Text::NSP::Measures::2D::Fisher::left;
use Text::NSP::Measures::2D::Fisher::twotailed;
use Statistics::Multtest qw(BH);
use SVG;
use SVG::Extension;
use SVG::DOM;
use SVG::XML;
use SVG::Element;
my $version="1.0 version";
my $program=$0;
use Getopt::Long;
my %opts;

my $usage=<<USAGE; #******* Instruction of this program *********#

	Author: Wanfei Liu
	Email: <liuwf\@big.ac.cn>
	Date: May 16, 2016
	Version: $version

	Introduction
	This program can assemble and/or annotate genome for new genome and known genome upgrade using sequence alignment file (SAM or BAM format), sequence variant file (VCF format or five coloum table (tab-delimited, including chromosome, position, id, reference allele and alternative allele)) or new genome sequence file (FASTA format) based on reference genome sequence file (FASTA format) and/or annotation file (TBL, GTF, GFF, GFF3 or BED format).
	
	Usage: perl $program -g genome_sequence(FASTA) -a annotation(TBL, GTF, GFF, GFF3 or BED) [-s SAM_file | -b BAM_file | -v VCF_file | -n new_genome_file] -o prefix_of_output_file.

	   -g: the absolute path of genome sequence file (FASTA format, required).
	   -a: the absolute path of genome annotation file (TBL, GTF, GFF, GFF3 or BED format, required for annotation). 
	   -s: the absolute path of SAM file (sorted according to coordinate, it is mututally exclusive with -b, -v and -n).
	   -b: the absolute path of BAM file (sorted according to coordinate, it is mututally exclusive with -s, -v and -n).
	   -v: the absolute path of sequence variant file (VCF format or five coloum table (tab-delimited, including chromosome, position, id, reference allele and alternative allele), it is mututally exclusive with -s, -b and -n).
	   -n: the absolute path of new genome sequence file (FASTA format, it is mututally exclusive with -s, -b and -v).
	   -m: the minIdentity for BLAT comparative between reference and new genome (default value is 90).
	   -q: quality threshold for reads (default value is 30).
	   -l: read length threshold (default value is 30).
	   -d: read depth threshold for sequence variant (default value is 3).
	   -c: allele depth threshold for sequence variant (default value is 3).
	   -p: allele proportion threshold for sequence variant (default value is 0.5).
	   -t: thread for parallel Blat
	  -fu: filter unpaired reads (yes or no, default value is no).
	  -fm: filter multiple mapping reads (yes or no, default value is yes)	
	   -o: the prefix of output file (required).

	Note: Genome sequence, annotation and sequence variant files should be uncompressed files. The name of annotation file should be ended with suffix *.tbl, *.gtf, *.gff, *.gff2, *.gff3 or *.bed. If you want run program multiple times with different parameters, please use "-o" option to assign different prefix name for output files. To filter false positive sequence variants, we use 3 as the default minimum reads depth and 3 as the default minimum alt allele reads coverage, because most of false positive variants come from the low depth and low alt coverage records according to our previous study. Although we filtered genome comparative result (obtained by BLAT) for reference and new genome sequence, user can manually check the *.psl file for higher accuracy if possible and only remain the most possible direction results (plus or minus). Moreover, if comparing with different assembly version or different strain of same species, we can use the default minIdentity parameter (90) for BLAT, however, if comparing with different species, we recommend using 50 as minIdentity by option -m (-m 50) for BLAT. To convert BAM file to SAM file, RGAAT need samtools.

USAGE

#Gather input
GetOptions(\%opts,"g:s","a:s","s:s","b:s","v:s","n:s","q:s","l:s","d:s","c:s","p:s","fu:s","fm:s","m:s","o:s","t:s");

#Verify input
if (!defined $opts{g} || (!defined $opts{s} and !defined $opts{b} and !defined $opts{v} and !defined $opts{n}) || !defined $opts{o}) {
	die "$usage\n";
}

my $genomefile=$opts{g};
my $annfile=$opts{a};
my $samfile=$opts{s};
my $bamfile=$opts{b};
my $vcffile=$opts{v};
my $genomenew=$opts{n};
my $prefix=$opts{o};
my $thread=$opts{t};
#Default parameters
my $quality = (defined $opts{q})?$opts{q}:30; #quality_threshold
my $read_len = (defined $opts{l})?$opts{l}:30; #read_length_threshold
my $depth_cutoff = (defined $opts{d})?$opts{d}:3;  #Default:3
my $alt_cutoff = (defined $opts{c})?$opts{c}:3; #Default:3
my $alt_proportion=(defined $opts{p})?$opts{p}:0.5; #Default minimal alt proportion is 0.1
my $filter_unpair=(defined $opts{fu})?$opts{fu}:"no"; #default value is no
my $filter_mult=(defined $opts{fm})?$opts{fm}:"yes"; #default value is yes
my $minidentity=(defined $opts{m})?$opts{m}:90; #default value is 90

#mismatch cutoff for read mapping result
my %miscut = (
	"0" => "0",
	"1" => "1",
	"2" => "3",
	"3" => "20",
	"4" => "47",
	"5" => "79",
	"6" => "114",
);
for (my $i=7;$i<500;$i++) {
	$miscut{$i}=$miscut{$i-1}+35;
}

my $cmdline = "$0";
foreach my $opt (keys %opts) {
	$cmdline.=" -$opt $opts{$opt}";
}

#path for program and depended files
my @scriptpath = split /\//, $0;
my $scriptname = pop @scriptpath;
my $scriptdir  = join '/', @scriptpath;

my (%change,%ref,%var,%pos,%seq,%ref_len,%new_len);

my $time=localtime;
print "Start time: $time.\n";
print "$cmdline\n";

if (exists $opts{s} or exists $opts{b}) {
	#identify sequence variants
	open (OUT,">$prefix.var")||die("fail to open $prefix.var.\n");
	print OUT "##fileformat=VCF like format\n";
	print OUT "##INFO=<ID=RDP,Number=1,Type=Integer,Description=\"All mapped reads depth after reads filter\">\n";
	print OUT "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"High quality reads depth; some reads and loci may have been filtered\">\n";
	print OUT "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency, for each ALT allele, in the same order as listed\">\n";
	print OUT "##INFO=<ID=BQ,Number=1,Type=Integer,Description=\"Locus base Quality\">\n";
	print OUT "##INFO=<ID=MP,Number=1,Type=Integer,Description=\"Multiple mapped reads depth\">\n";
	print OUT "##INFO=<ID=IC,Number=1,Type=Integer,Description=\"Insertion reads count\">\n";
	print OUT "##INFO=<ID=DC,Number=1,Type=Integer,Description=\"Deletion reads count\">\n";
	print OUT "##INFO=<ID=PC,Number=1,Type=Integer,Description=\"Plus strand reads count\">\n";
	print OUT "##INFO=<ID=MC,Number=1,Type=Integer,Description=\"Minus strand reads count\">\n";
	print OUT "##INFO=<ID=SC,Number=1,Type=Integer,Description=\"Startpoint count\">\n";
	print OUT "##INFO=<ID=Ref,Number=1,Type=String,Description=\"Reference allele\">\n";
	print OUT "##INFO=<ID=Alt,Number=1,Type=String,Description=\"Alternative allele\">\n";
	&sam2variant;
	close OUT;
	%ref=();%var=();
	#local assembly
	##update genome and annotation
	my $genomeupdate="$prefix.genome_update";
	my $poschange="$prefix.pos_change";
	my $annoupdate="$prefix.anno_update";
	open (GENOME,">$genomeupdate")||die("fail to open $genomeupdate.\n");
	open (POS,">$poschange")||die("fail to open $poschange.\n");
	print POS "#Chr\tOldStart\tOldEnd\tNewStart\tNewEnd\tStrand\n";
	open (ANNO, ">$annoupdate")||die("fail to open $annoupdate.\n") if (exists $opts{a});
	
	$time=localtime;
	print "Read genome file: $time.\n";
	open (IN,"<$genomefile")||die("fail to open $genomefile.\n");
	my $id=undef;
	while (<IN>) {
		chomp;
		if (/^>(\S+)/) {
			$id=$1;
			next;
		}
		$ref{$id}.=$_;
	}
	close IN;

	$time=localtime;
	print "Update genome: $time.\n";
	foreach $id (keys %ref) {
		$ref_len{$id}=length($ref{$id});
		&genome_update($id,$ref{$id});
	}
	close GENOME;
	close POS;
	close ANNO if (exists $opts{a});

	$time=localtime;
	print "Genome comparison: $time.\n";
	&genome_compare_figure;
}elsif (exists $opts{v}) {
	$time=localtime;
	print "Reads sequence variant file: $time.\n";
	#identify sequence variants
	&variant;

	#update genome and annotation
	my $genomeupdate="$prefix.genome_update";
	my $poschange="$prefix.pos_change";
	my $annoupdate="$prefix.anno_update";

	open (GENOME,">$genomeupdate")||die("fail to open $genomeupdate.\n");
	open (POS,">$poschange")||die("fail to open $poschange.\n");
	print POS "#Chr\tOldStart\tOldEnd\tNewStart\tNewEnd\tStrand\n";
	open (ANNO, ">$annoupdate")||die("fail to open $annoupdate.\n") if (exists $opts{a});

	$time=localtime;
	print "Read genome file: $time.\n";
	open (IN,"<$genomefile")||die("fail to open $genomefile.\n");
	my $id=undef;
	while (<IN>) {
		chomp;
		if (/^>(\S+)/) {
			$id=$1;
			next;
		}
		$ref{$id}.=$_;
	}
	close IN;

	$time=localtime;
	print "Update genome: $time.\n";
	foreach $id (keys %ref) {
		$ref_len{$id}=length($ref{$id});
		&genome_update($id,$ref{$id});
	}
	close GENOME;
	close POS;
	close ANNO if (exists $opts{a});

	$time=localtime;
	print "Genome comparison: $time.\n";
	&genome_compare_figure;
}elsif (exists $opts{n}) {
	#genome comparative
	# Assembly and strain: nohup blat database query -minIdentity=90 -noHead out.psl
	# Species: nohup blat database query -minIdentity=50 -noHead out.psl

	$time=localtime;
	print "Genome comparison by BLAT: $time.\n";
	system "./blat_smp.pl -t $genomefile -q $genomenew -o $prefix.psl -p $thread -- -minIdentity=$minidentity -minScore=100 -noHead";
	#blat result filter

	$time=localtime;
	print "BLAT result filter: $time.\n";
	&blat_filter("$prefix.psl","$prefix.filter.psl");
	#identify sequence variants and create coordinate transform 

	$time=localtime;
	print "Identify sequence variant: $time.\n";
	&variant_blat("$prefix.filter.psl");
	#annotation_transfer

	$time=localtime;
	print "Annotation transfer: $time.\n";
	&annotation_transfer;

	$time=localtime;
	print "Genome comparison: $time.\n";
	&genome_compare_figure;
}

#local assembly

$time=localtime;
print "End time: $time.\n";

sub sam2variant {
	$time=localtime;
	print "Start identify sequence variant: $time.\n";
	my $freq_cut=$alt_proportion;
	my $cov_var=$alt_cutoff;
	my $cov_read=$depth_cutoff;

	my $read_num=0;
	my $read_filtered=0;
	my $filter_unmap=0;
	my $filter_pcr=0;
	my $filter_quality=0;
	my $filter_mult=0;
	my $filter_short=0;
	my $filter_qual=0;
	my $startpos=0;
	my $chrseq=undef;
	my $id=undef;
	if (exists $opts{s}) {
		open (STDIN,"<$opts{s}")||die("fail to open $opts{s}.\n");
	}elsif (exists $opts{b}) {
		open STDIN,"samtools view $opts{b} |";
	}
	while (<STDIN>) {
		chomp;
		#SRR027895.5259200_1     99      chr1    15443   255     76M     =       15623   256     CTCCAGAGGCCTCAGGTCCAGTCTCTAAAAATATCTCAGGAAGCTGCAGTGGCTGACCATTGCCTTGGACCGCTCT    BCBCCBBABCBBBBBB=BBAB?BBBBBAABBBBBBBBBBA@BBBABBBB?@AAA>?<702>7:;1:+670:9;<==    NM:i:1  NH:i:1
		## Bit    Description                                                Comment                                Value
		## 0x1    template having multiple segments in sequencing            0: single-end 1: paired end            value: 2^^0 (  1)
		## 0x2    each segment properly aligned according to the aligner     true only for paired-end alignments    value: 2^^1 (  2)
		## 0x4    segment unmapped                                           ---                                           ---
		## 0x8    next segment in the template unmapped                      ---                                           ---
		## 0x10   SEQ being reverse complemented                                                                    value: 2^^4 ( 16)
		## 0x20   SEQ of the next segment in the template being reversed                                            value: 2^^5 ( 32)
		## 0x40   the first segment in the template                          read 1                                 value: 2^^6 ( 64)
		## 0x80   the last segment in the template                           read 2                                 value: 2^^7 (128)
		## 0x100  secondary alignment                                        ---                                           ---
		## 0x200  not passing quality controls                               ---                                           ---
		## 0x400  PCR or optical duplicate                                   ---                                           ---
		
		next if(/^\@/);
		$read_num++;
		if ($read_num%1000000==0) {
			print "Processed $read_num reads.\n";
		}
		my @list=split /\t/,$_;

		#read chromosome
		if (!defined $id or $id ne $list[2]) {
			&vcf_print($id) if (defined $id);
			$chrseq=&genomeseq($genomefile,$list[2]);
			$id=$list[2];%var=();%ref=();
		}

		#filter reads by average base quality
		my $qual=0;
		if ($list[10] ne "*") {
			foreach my $basequal (split //,$list[10]) {
				$qual+=ord($basequal)-33;
			}
			$qual=$qual/(length($list[10]));
			if ($qual<$quality) {
				$filter_qual++;
				next;
			}
		}

		#filter read according to the mapping quality marked in flag
		my $bin=unpack("B32",pack("N",$list[1]));
		$bin=sprintf "%011d",$bin;
		my @bin=split //,$bin;
		if ($opts{fu} ne "no") {
			if (($bin[7] eq "1") or ($bin[8] eq "1")) {
				$filter_unmap++;
				next;
			}
		}else {
			if ($bin[8] eq "1") {
				$filter_unmap++;
				next;
			}
		}
		if ($bin[0] eq "1") {
			$filter_pcr++;
			next;
		}
		if ($bin[1] eq "1") {
			$filter_quality++;
			next;
		}
		if ($opts{fm} ne "no") {
			if ($bin[2] eq "1") {
				$filter_mult++;
				next;
			}
		}

		#filter reads by mismatch, insertion and deletion
		my $nm=0;
		my @md=();
		if (/NM:i:(\d+)/) {
			$nm=$1;
		}
		if (/MD:Z:([\S+])/) {
			@md=($1=~m/(\D+)/g);
		}
		if ($miscut{scalar(@md)}>length($list[9])) {
			$filter_quality++;
			next;
		}elsif ($miscut{$nm}>length($list[9])) {
			$filter_quality++;
			next;
		}
		
		#clip read ends with short match regions
		$list[5]=~s/\d+H+//g;
		my @record=($list[5]=~m/(\d+\D+)/g);
		my @num=();
		my @label=();
		my $SM=undef;
		my $EM=undef;
		my $matchlen=0;
		foreach my $record (@record) {
			if ($record=~/(\d+)(\D+)/) {
				push (@num,$1);
				push (@label,$2);
				$matchlen+=$1 if ($2 eq "M" or $2 eq "=" or $2 eq "I");
				$SM=$1 if (!defined $SM and ($2 eq "M" or $2 eq "="));
				$EM=$1 if ($2 eq "M" or $2 eq "=");
			}
		}

		#filter read length
		if ($matchlen<$read_len) {
			$filter_short++;
			next;
		}

		while ($SM<8) {
			my ($SMnew,$clipv,$clipr,$label,$num)=&headclip (\@label,\@num);
			$SM=$SMnew;
			@label=@{$label};
			@num=@{$num};
			if (@label<1) {
				last;
			}else {
				$list[3]=$list[3]+$clipr;
				$list[9]=substr($list[9],$clipv,length($list[9]));
				$list[10]=substr($list[10],$clipv,length($list[10]));
			}
		}
		if (@label<1) {
			$filter_short++;
			next;
		}
		while ($EM<8) {
			my ($EMnew,$clipv,$label,$num)=&tailclip (\@label,\@num);
			$EM=$EMnew;
			@label=@{$label};
			@num=@{$num};
			if (@label<1) {
				last;
			}else {
				$list[9]=substr($list[9],0,length($list[9])-$clipv);
				$list[10]=substr($list[10],0,length($list[10])-$clipv);
			}
		}
		if (@label<1) {
			$filter_short++;
			next;
		}

		#process reads
		$read_filtered++;
		my $clip=0;
		my $startr=0;
		my $startv=0;
		for (my $i=0;$i<@label;$i++) {
			if ($label[$i] eq "S") {
				$clip=$num[$i];
			}elsif ($label[$i] eq "M") {
				for (my $j=0;$j<$num[$i];$j++) {
					$var{$list[3]+$startr}{"all"}++;
					my $qua=substr($list[10],$clip+$startv,1);
					$qua=ord($qua)-33;
					if ($qua>=0.85*$qual) {
						my $var=substr($list[9],$clip+$startv,1);
						my $ref=substr($chrseq,$list[3]+$startr-1,1);
						if ($ref=~/[acgtnN]+/ or $var=~/[nN]+/) {
						}else {
							$var=uc($var);$ref=uc($ref);
							$var{$list[3]+$startr}{"total"}++;
							$var{$list[3]+$startr}{"qual"}+=$qua;
							$var{$list[3]+$startr}{"start"}++ if ($list[3] ne $startpos);
							if ($bin[2] eq "1") {
								$var{$list[3]+$startr}{"mult"}++;
							}
							if ($bin[6] eq "1") {
								$var{$list[3]+$startr}{"minus"}++;
							}else {
								$var{$list[3]+$startr}{"plus"}++;
							}
							if ($ref eq $var) {
								$var{$list[3]+$startr}{"ref"}++;
							}else {
								$var{$list[3]+$startr}{"$ref\t$var"}++;
							}
						}
					}
					$startr++;
					$startv++;
				}
			}elsif ($label[$i] eq "N") {
				$startr+=$startr+$num[$i];
			}elsif ($label[$i] eq "I") {
				my $qua_str=substr($list[10],$clip+$startv,$num[$i]);
				my $qua=0;
				foreach my $basequal (split //,$qua_str) {
					$qua+=ord($basequal)-33;
				}
				$qua=$qua/($num[$i]);
				if ($qua>=0.85*$qual) {
					my $var=substr($list[9],$clip+$startv-1,$num[$i]+1);
					my $ref=substr($chrseq,$list[3]+$startr-2,1);
					if ($ref=~/[acgtnN]+/ or $var=~/[nN]+/) {
					}else {
						$var=uc($var);$ref=uc($ref);
						$var{$list[3]+$startr-1}{"$ref\t$var"}++;
						$var{$list[3]+$startr-1}{"insert"}++;
					}
				}
				$startv+=$num[$i];
			}elsif ($label[$i] eq "D") {
				my $qua_str=substr($list[10],$clip+$startv-1,2);
				my $qua=0;
				foreach my $basequal (split //,$qua_str) {
					$qua+=ord($basequal)-33;
				}
				$qua=$qua/2;
				if ($qua>=0.85*$qual) {
					my $var=substr($list[9],$clip+$startv-1,1);
					my $ref=substr($chrseq,$list[3]+$startr-2,$num[$i]+1);
					if ($ref=~/[acgtnN]+/ or $var=~/[nN]+/) {
					}else {
						$var=uc($var);$ref=uc($ref);
						$var{$list[3]+$startr-1}{"$ref\t$var"}++;
						$var{$list[3]+$startr-1}{"delete"}++;
					}
				}
				$startr+=$num[$i];
			}
		}
		#start position
		if ($list[3] ne $startpos) {
			$startpos=$list[3];
		}
	}
	close STDIN;
	#last chromosome
	&vcf_print($id);

	print "#Total reads is $read_num.\n#Total filtered reads is $read_filtered\n";
	print "#The filter reads are as following:\n";
	print "Unmap\t$filter_unmap\n";
	print "PCR_duplicate\t$filter_pcr\n";
	print "Quality_control\t$filter_quality\n";
	print "Multiple\t$filter_mult\n";
	print "Short_anchor\t$filter_short\n";
	print "Low_quality\t$filter_qual\n";

	$time=localtime;
	print "Sequence variant identified: $time.\n";
}

sub vcf_print {
	my $id=shift;

	my $indelf=0;#the frequency of before indel
	my $indelp=0;#the position of before indel
	my $indell=0;#the length of before indel

	foreach my $pos (sort {$a <=> $b} keys %var) {
		my $all=$var{$pos}{"all"};
		my $total=$var{$pos}{"total"};
		my $qual=0;
		if (defined $var{$pos}{"qual"}) {
			$qual=(sprintf "%.0f",$var{$pos}{"qual"}/$total) if ($total>0);
		}
		my $refcount=0;
		if (defined $var{$pos}{"ref"}) {
			$refcount=$var{$pos}{"ref"};
		}
		my $multcount=0;
		if (defined $var{$pos}{"mult"}) {
			$multcount=$var{$pos}{"mult"};
		}
		my $insertcount=0;
		if (defined $var{$pos}{"insert"}) {
			$insertcount=$var{$pos}{"insert"};
		}
		my $deletecount=0;
		if (defined $var{$pos}{"delete"}) {
			$deletecount=$var{$pos}{"delete"};
		}
		my $pluscount=0;
		if (defined $var{$pos}{"plus"}) {
			$pluscount=$var{$pos}{"plus"};
		}
		my $minuscount=0;
		if (defined $var{$pos}{"minus"}) {
			$minuscount=$var{$pos}{"minus"};
		}
		my $startcount=0;
		if (exists $var{$pos}{"start"}) {
			$startcount=$var{$pos}{"start"};
		}
		if ($total<$depth_cutoff) {
			next;
		}
		my $flags=0;
		my $flagi=0;
		foreach my $var (sort { $var{$pos}{$b}<=>$var{$pos}{$a} } keys %{$var{$pos}}) {
			next if ($var eq "all");
			next if ($var eq "total");
			next if ($var eq "qual");
			next if ($var eq "ref");
			next if ($var eq "mult");
			next if ($var eq "insert");
			next if ($var eq "delete");
			next if ($var eq "plus");
			next if ($var eq "minus");
			next if ($var eq "start");
			my ($ref,$read)=split /\t/,$var;
			next if ($var{$pos}{$var}<$alt_cutoff);
			my $varcount=$var{$pos}{$var};
			my $AF=(sprintf "%.2f",$varcount/$total);
			next if ($AF<$alt_proportion);

			if (length($ref)==length($read) and $flags==0 and $flagi==0) {
				if ($pos<=$indelp+$indell-1) {
					if ($indelf<$AF) {
						if ($multcount>0) {
							if (($varcount-$multcount)/$total>=$alt_proportion) {
								if (($pluscount>0 and $minuscount>0) or $startcount>=2) {
									printf OUT "$id\t$pos\t\.\t$ref\t$read\t\.\t\.\tRDP=$all;DP=$total;AF=$AF;BQ=$qual;MP=$multcount;IC=$insertcount;DC=$deletecount;PC=$pluscount;MC=$minuscount;SC=$startcount;Ref=$refcount;Alt=$varcount\n";
									$change{$id}{$pos}="$ref\t$read";
								}
							}
						}else {
							if (($pluscount>0 and $minuscount>0) or $startcount>=2) {
								printf OUT "$id\t$pos\t\.\t$ref\t$read\t\.\t\.\tRDP=$all;DP=$total;AF=$AF;BQ=$qual;MP=$multcount;IC=$insertcount;DC=$deletecount;PC=$pluscount;MC=$minuscount;SC=$startcount;Ref=$refcount;Alt=$varcount\n";
								$change{$id}{$pos}="$ref\t$read";
							}
						}
					}
				}else {
					if ($multcount>0) {
						if (($varcount-$multcount)/$total>=$alt_proportion) {
							if (($pluscount>0 and $minuscount>0) or $startcount>=2) {
								printf OUT "$id\t$pos\t\.\t$ref\t$read\t\.\t\.\tRDP=$all;DP=$total;AF=$AF;BQ=$qual;MP=$multcount;IC=$insertcount;DC=$deletecount;PC=$pluscount;MC=$minuscount;SC=$startcount;Ref=$refcount;Alt=$varcount\n";
								$change{$id}{$pos}="$ref\t$read";
							}
						}
					}else {
						if (($pluscount>0 and $minuscount>0) or $startcount>=2) {
							printf OUT "$id\t$pos\t\.\t$ref\t$read\t\.\t\.\tRDP=$all;DP=$total;AF=$AF;BQ=$qual;MP=$multcount;IC=$insertcount;DC=$deletecount;PC=$pluscount;MC=$minuscount;SC=$startcount;Ref=$refcount;Alt=$varcount\n";
							$change{$id}{$pos}="$ref\t$read";
						}
					}
				}
				$flags=1;
			}

			if (length($ref)!=length($read) and $flagi==0) {
				my $indelcount=0;#ref count of indel
				if (length($ref)<length($read)) {
					$indelcount=$total-$insertcount if ($total-$insertcount>0);
				}else {
					$indelcount=$total-$deletecount if ($total-$deletecount>0);
				}
				my $flag=0;
				if ($multcount>0) {
					if (($varcount-$multcount)/$total>=$alt_proportion) {
						if (($pluscount>0 and $minuscount>0) or $startcount>=2) {
							printf OUT "$id\t$pos\t\.\t$ref\t$read\t\.\t\.\tRDP=$all;DP=$total;AF=$AF;BQ=$qual;MP=$multcount;IC=$insertcount;DC=$deletecount;PC=$pluscount;MC=$minuscount;SC=$startcount;Ref=$indelcount;Alt=$varcount\n";
							$change{$id}{$pos}="$ref\t$read";
							$flag=1;
						}
					}
				}else {
					if (($pluscount>0 and $minuscount>0) or $startcount>=2) {
						printf OUT "$id\t$pos\t\.\t$ref\t$read\t\.\t\.\tRDP=$all;DP=$total;AF=$AF;BQ=$qual;MP=$multcount;IC=$insertcount;DC=$deletecount;PC=$pluscount;MC=$minuscount;SC=$startcount;Ref=$indelcount;Alt=$varcount\n";
						$change{$id}{$pos}="$ref\t$read";
						$flag=1;
					}
				}
				$flagi=1;
				if ($flag=1 and length($read)<length($ref)) {
					$indelf=$var{$pos}{$var}/$total;
					$indell=length($ref);
					$indelp=$pos;
				}
			}
			last if ($flags==1 and $flagi==1);
		}
	}
	return;
}

sub headclip {
	my $label=shift;
	my @label=@{$label};
	my $num=shift;
	my @num=@{$num};

	my $clipv=0;
	my $clipr=0;
	my @numnew=();
	my @labelnew=();
	my $SM=undef;
	my $M=0;
	for (my $i=0;$i<@label;$i++) {
		if ($label[$i] eq "S") {
			if ($M<1) {
				$clipv+=$num[$i];
			}else {
				push (@numnew,$num[$i]);
				push (@labelnew,$label[$i]);
			}
		}elsif ($label[$i] eq "M") {
			if ($M<1) {
				$clipv+=$num[$i];
				$clipr+=$num[$i];
				$M++;
			}else {
				$SM=$num[$i] if (!defined $SM);
				push (@numnew,$num[$i]);
				push (@labelnew,$label[$i]);
				$M++;
			}
		}elsif ($label[$i] eq "N") {
			if ($M<1) {
				$clipr+=$num[$i];
			}else {
				push (@numnew,$num[$i]);
				push (@labelnew,$label[$i]);
			}
		}elsif ($label[$i] eq "I") {
			if ($M<1) {
				$clipv+=$num[$i];
			}else {
				push (@numnew,$num[$i]);
				push (@labelnew,$label[$i]);
			}
		}elsif ($label[$i] eq "D") {
			if ($M<1) {
				$clipr+=$num[$i];
			}else {
				push (@numnew,$num[$i]);
				push (@labelnew,$label[$i]);
			}
		}
	}
	@label=@labelnew;
	@num=@numnew;
	return ($SM,$clipv,$clipr,\@label,\@num);
}

sub tailclip {
	my $label=shift;
	my @label=@{$label};
	my $num=shift;
	my @num=@{$num};

	my $clipv=0;
	my @numnew=();
	my @labelnew=();
	my $EM=undef;
	my $M=0;

	for (my $i=$#label;$i>=0;$i--) {
		if ($label[$i] eq "S") {
			if ($M<1) {
				$clipv+=$num[$i];
			}else {
				unshift (@numnew,$num[$i]);
				unshift (@labelnew,$label[$i]);
			}
		}elsif ($label[$i] eq "M") {
			if ($M<1) {
				$clipv+=$num[$i];
				$M++;
			}else {
				$EM=$num[$i] if (!defined $EM);
				unshift (@numnew,$num[$i]);
				unshift (@labelnew,$label[$i]);
				$M++;
			}
		}elsif ($label[$i] eq "N") {
			if ($M<1) {
			}else {
				unshift (@numnew,$num[$i]);
				unshift (@labelnew,$label[$i]);
			}
		}elsif ($label[$i] eq "I") {
			if ($M<1) {
				$clipv+=$num[$i];
			}else {
				unshift (@numnew,$num[$i]);
				unshift (@labelnew,$label[$i]);
			}
		}elsif ($label[$i] eq "D") {
			if ($M<1) {
			}else {
				unshift (@numnew,$num[$i]);
				unshift (@labelnew,$label[$i]);
			}
		}
	}
	@label=@labelnew;
	@num=@numnew;
	return ($EM,$clipv,\@label,\@num);
}

sub genomeseq {
	my $infile=shift;
	my $chr=shift;
	open (IN,"<$infile")||die("fail to open $infile.\n");
	my $id=undef;
	my $seq=undef;
	while(<IN>){
		next if(/^\#/); #ignore header
		chomp;
		if (/>(\S+)/) {
			$id=$1;
			last if (defined $seq);
			next;
		}
		if ($id ne $chr) {
			next;
		}else {
			$seq.=$_;
		}
	}
	close IN;
	return $seq;
}

sub variant {
	my $pos=undef;
	my $len=undef;
	open (IN,"<$opts{v}")||die("fail to open $opts{v}.\n");
	while (<IN>) {
		chomp;
		#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SRR1063404_cp.sorted.bam
		#mt      30088   .       C       T       .       .       DP=4904;AF=0.17;AC=821
		#mt      30101   .       A       AT      .       .       INDEL;IDV=3665;IMF=0.86;DP=4255
		next if (/^\#/);
		$_=~s/\,?\<NON\_REF\>//;
		my @list=split /\t/,$_;
		next if ($list[7]=~/^END\=\d+/);
		next if ($list[4] eq ".");
		my %attr=();
		my @ALT=split /\,/,$list[4];
		for (my $j=0;$j<@ALT;$j++) {
			$attr{"ALT$j"}=$ALT[$j];
		}
		if (length($list[3])!=length($attr{"ALT0"})) {
			if (!defined $pos) {
				$pos=$list[1];
				$len=length($list[3]);
			}else {
				if ($list[1]<$pos+$len-1) {
					$pos=$list[1];
					$len=length($list[3]);
					next;
				}else {
					$pos=$list[1];
					$len=length($list[3]);
				}
			}
			my $clip1=0;
			my $clip2=0;
			my @ref=split //,$list[3];
			my @var=split //,$list[4];
			for (my $i=0;$i<@ref-1 and $i<@var-1;$i++) {
				if ($ref[$#ref-$i] eq $var[$#var-$i]) {
					$clip2++;
				}else {
					last;
				}
			}
			for (my $i=0;$i<@ref-$clip2-1 and $i<@var-$clip2-1;$i++) {
				if ($ref[$i] eq $var[$i]) {
					$clip1++;
				}else {
					last;
				}
			}
			if ($clip1>0) {
				$list[3]=substr($list[3],$clip1-1,length($list[3])-($clip1-1)-$clip2);
				$list[4]=substr($list[4],$clip1-1,length($list[4])-($clip1-1)-$clip2);
				$list[1]=$list[1]+$clip1-1;
			}else {
				$list[3]=substr($list[3],0,length($list[3])-$clip2);
				$list[4]=substr($list[4],0,length($list[4])-$clip2);
			}
		}
		if (defined $list[7]) {
			my @attributes = split /;/, $list[7];
			foreach my $attr ( @attributes) {
				next unless $attr =~ /^(\S+)\=(\S+)$/;
				my $c_type  = $1;
				my $c_value = $2;
				$attr{$c_type} = $c_value;
			}
			my $ref=0;
			my $alt=0;
			if ($list[7]=~/DP4\=(\d+)\,(\d+)\,(\d+)\,(\d+)/) {
				$ref=$1+$2;
				$alt=$3+$4;
			}elsif ($list[8]=~/AD/) {
				my @attributes1 = split /:/, $list[8];
				my @attributes2 = split /:/, $list[9];
				for (my $j=0;$j<@attributes1;$j++) {
					my $c_type  = $attributes1[$j];
					my $c_value = $attributes2[$j];
					$attr{$c_type} = $c_value;
					my @c_value=split /\,/,$c_value;
					for (my $k=0;$k<@c_value;$k++) {
						$attr{"$c_type$k"}=$c_value[$k];
					}
				}
				$ref=$attr{"AD0"};
				$alt=$attr{"AD1"};
			#Ref=1744;Alt=3967
			}elsif ($list[7]=~/Ref\=([^\;]+)\;Alt\=(\d+)/) {
				$ref=$1 if ($1>0);
				$alt=$2;
			}
			next if ($alt+$ref<$depth_cutoff);
			my $AF=$alt/($alt+$ref);
			next if ($alt/($ref+$alt)<$alt_proportion);
			$change{$list[0]}{$list[1]}="$list[3]\t$list[4]";
		}else {
			$change{$list[0]}{$list[1]}="$list[3]\t$list[4]";
		}
	}
	close IN;
	return;
}

sub genome_update {	
	my $id=shift;
	my $seq=shift;

	my $old_start=undef;
	my $old_end=0;
	my $new_start=undef;
	my $new_end=0;
	my $str=undef;
	my $num=0;
	my $total=0;
	%pos=();
	for (my $i=0;$i<length($seq);$i++) {
		$old_start=$i+1 if (!defined $old_start);
		$old_end=$i+1;
		$new_start=$total+1 if (!defined $new_start);
		$new_end=$total+1;
		if (!defined $change{$id}{$i+1}) {
			$str.=substr($seq,$i,1);
			$pos{$i+1}=$total+1;
			$total++;
		}else {
			print POS "$id\t$old_start\t$old_end\t$new_start\t$new_end\tNewGenome\t+\n";
			$old_start=undef;
			$new_start=undef;
			my ($ref,$change)=split /\t/,$change{$id}{$i+1};
			$str.=$change;
			$pos{$i+1}=$total+1;			
			my $num=length($change)-length($ref);
			if ($num>=0) {
				for (my $j=0;$j<$num;$j++) {
					$total++;
				}
			}else {
				for (my $j=0;$j<abs($num);$j++) {
					$pos{$i+1+$j+1}=$total+1;
				}
			}
			$total++;
			$i=$i+length($ref)-1;
		}
	}
	print POS "$id\t$old_start\t$old_end\t$new_start\t$new_end\tNewGenome\t+\n";
	print GENOME ">$id\n";
	$new_len{$id}=length($str);
	for ( my $pos = 0 ; $pos < length($str) ; $pos += 60 ) {
		print GENOME (substr($str, $pos, 60))."\n";
	}
	#print annotation
	next if (!exists $opts{a});
	if ($annfile=~/\.tbl$/) {
		&tblupdate($annfile,\%pos,$id);
	}elsif ($annfile=~/\.gtf$/ or $annfile=~/\.gff$/ or $annfile=~/\.gff2$/ or $annfile=~/\.gff3$/) {
		&gffupdate($annfile,\%pos,$id);
	}elsif ($annfile=~/\.bed$/) {
		&bedupdate($annfile,\%pos,$id);
	}
	return;
}

sub tblupdate {
	my $infile=shift;
	my $pos=shift;
	my %pos=%{$pos};
	my $chr=shift;

	my $id;
	my $flag=0;
	open (IN,"<$infile")||die("fail to open $infile.\n");
	while (<IN>) {
		chomp;
		if (/^>\S+\s+(\S+)$/) {
			$id=$1;
			if ($id eq $chr) {
				$flag=1;
			}else {
				$flag=0;
			}
			print ANNO "$_\n" if ($flag==1);
		}elsif (/^<?(\d+)\t(\d+)\t(\S+)$/) {
			#122026460	125184587	centromere
			#11874	14409	gene
			my $start="?";
			$start=$pos{$1} if (exists $pos{$1});
			my $end="?";
			$end=$pos{$2} if (exists $pos{$2});
			print ANNO "$start\t$end\t$3\n" if ($flag==1);
		}elsif (/^\t\t\t\S+\t.+$/) {
			#			product	DEAD/H (Asp-Glu-Ala-Asp/His) box helicase 11 like 1
			#			note	Linear centromere model derived predominantly from reads generated in PMID: 17803354. This region does not represent an actual centromere sequence, as long-range ordering of repeats and unmapped WGS contigs is not provided by the model. For details of model production, see http://arxiv.org/abs/1307.0035.
			#			gene	DDX11L1
			print ANNO "$_\n" if ($flag==1);
		}elsif (/^<?(\d+)\t(\d+)$/) {
			my $start="?";
			$start=$pos{$1} if (exists $pos{$1});
			my $end="?";
			$end=$pos{$2} if (exists $pos{$2});
			print ANNO "$start\t$end\n" if ($flag==1);
		}else {
			print ANNO "$_\n" if ($flag==1);
		}
	}
	close IN;
	return;
}

sub gffupdate {
	my $infile=shift;
	my $pos=shift;
	my %pos=%{$pos};
	my $chr=shift;

	open (IN,"<$infile")||die("fail to open $infile.\n");
	while (<IN>) {
		chomp;
		#chr1	Cufflinks	exon	4886744	4886831	.	+	.	gene_id "XLOC_000005"; transcript_id "TCONS_00000014"; exon_number "4"; oId "CUFF.13.9"; tss_id "TSS8";
		next if (/^\#/);
		my @list=split/\t/,$_;
		next if ($list[0] ne $chr);
		my $start="?";
		$start=$pos{$list[3]} if (exists $pos{$list[3]});
		my $end="?";
		$end=$pos{$list[4]} if (exists $pos{$list[4]});
		my @listnew=(@list[0..2],$start,$end,@list[5..$#list]);
		my $listnew=join "\t",@listnew;
		print ANNO "$listnew\n";
	}
	close(IN);
	return;
}

sub bedupdate {
	my $infile=shift;
	my $pos=shift;
	my %pos=%{$pos};
	my $chr=shift;

	open (IN,"<$infile")||die("fail to open $infile.\n");
	while (<IN>) {
		chomp;
		#chr22 1000 5000 cloneA 960 + 1000 5000 0 2 567,488, 0,3512
		next if (/^\#/ or /^track/);
		my @list=split/\s+/,$_;
		next if ($list[0] ne $chr);
		my $start="?";
		$start=$pos{$list[1]+1}-1 if (exists $pos{$list[1]+1});
		my $end="?";
		$end=$pos{$list[2]} if (exists $pos{$list[2]});
		my @listnew=();
		if (@list>=10) {
			my $len;
			my $sta;
			my @len=split /\,/,$list[10];
			my @start=split /\,/,$list[11];
			for (my $i=0;$i<@start;$i++) {
				my $len_new="?";
				my $sta_new="?";
				$sta_new=$pos{$list[1]+$start[$i]+1}-$pos{$list[1]+1} if (exists $pos{$list[1]+1} and exists $pos{$list[1]+$start[$i]+1});
				$len_new=$pos{$list[1]+$start[$i]+$len[$i]}-$pos{$list[1]+$start[$i]+1}+1 if (exists $pos{$list[1]+$start[$i]+$len[$i]} and exists $pos{$list[1]+$start[$i]+1});
				$sta.=$sta_new.",";
				$len.=$len_new.",";
			}
			@listnew=($list[0],$start,$end,@list[3..5],$start,$end,@list[8..9],$len,$sta);
		}else {
			if (defined $list[3]) {
				@listnew=($list[0],$start,$end,@list[3..$#list]);
			}else {
				@listnew=($list[0],$start,$end);
			}
		}
		my $listnew=join "\t",@listnew;
		print ANNO "$listnew\n";
	}
	close(IN);
	return;
}

sub blat_filter {
	my $blatfile=shift;
	my $outfile=shift;

	my @list=();
	my $plusmatch=0;
	my $minusmatch=0;
	my $plusmis=0;
	my $minusmis=0;
	my $plusqueryinsert=0;
	my $minusqueryinsert=0;
	my $plustargetinsert=0;
	my $minustargetinsert=0;

	open (IN,"<$blatfile")||die("fail to open $blatfile.\n");
	while (<IN>) {
		chomp;
		#3640	232	0	0	21	240	18	153	+	gi|270610521|gb|GU183365.1|	37320	0	4112	S2846	19749	11158	15183	29	652,720,679,431,276,155,62,32,175,16,8,6,24,18,14,5,8,44,37,106,4,58,157,56,12,51,12,44,10,	0,661,1386,2095,2529,2806,2964,3052,3084,3259,3275,3288,3294,3325,3385,3399,3412,3432,3476,3514,3621,3629,3689,3846,3930,3950,4025,4040,4102,	11158,11810,12533,13212,13645,13921,14081,14144,14182,14360,14380,14388,14397,14421,14439,14456,14464,14510,14555,14599,14705,14714,14772,14952,15039,15065,15116,15129,15173,
		my @temp=split /\t/,$_;
		if ($temp[8] eq "+") {
			$plusmatch+=$temp[0];
			$plusmis+=$temp[1];
			$plusqueryinsert+=$temp[5];
			$plustargetinsert+=$temp[7];
		}else {
			$minusmatch+=$temp[0];
			$minusmis+=$temp[1];
			$minusqueryinsert+=$temp[5];
			$minustargetinsert+=$temp[7];
		}
		push @list,\@temp;
	}
	close IN;
	@list=sort {$a->[15] <=> $b->[15]} @list;
	@list=sort {$a->[13] cmp $b->[13]} @list;

	my $strand="+";
	if ( ($plusmatch+$plusmis)/($plusmatch+$plusmis+$plusqueryinsert+$plustargetinsert)<($minusmatch+$minusmis)/($minusmatch+$minusmis+$minusqueryinsert+$minustargetinsert) ) {
		$strand="-";
	}

	#remove bad alignment strand
	my @listnew=();
	for (my $i=0;$i<@list;$i++) {
		my @temp=@{$list[$i]};
		if ($temp[8] eq $strand) {
			push (@listnew,\@temp);
		}
	}
	@list=@listnew;
	@listnew=();

	open (OUT,">$outfile")||die("fail to open $outfile.\n");
	my @before=@{$list[0]};
	my $match=$before[0];
	my $mis=$before[1];
	my $queryinsert=$before[5];
	my $targetinsert=$before[7];
	my $start=$before[15];
	my $end=$before[16];
	for (my $i=1;$i<@list;$i++) {
		my @temp=@{$list[$i]};
		if ($temp[15]>$end) {
			print OUT join("\t",@before)."\n";
			$match=$temp[0];
			$mis=$temp[1];
			$queryinsert=$temp[5];
			$targetinsert=$temp[7];
			$start=$temp[15];
			$end=$temp[16];
			@before=@temp;
		}elsif ($start<=$temp[15] and $end>=$temp[16]) {
			if ( ($match+$mis)/($match+$mis+$queryinsert+$targetinsert)<($temp[0]+$temp[1])/($temp[0]+$temp[1]+$temp[5]+$temp[7]) ) {
				$match=$temp[0];
				$mis=$temp[1];
				$queryinsert=$temp[5];
				$targetinsert=$temp[7];
				$start=$temp[15];
				$end=$temp[16];
				@before=@temp;
			}
		}elsif ($end>=$temp[15] and $end<=$temp[16]) {
			if ( ($end-$temp[15])/($end-$start)<0.2 and ($end-$temp[15])/($temp[16]-$temp[15])<0.2 ) {
				print OUT join("\t",@before)."\n";
				$match=$temp[0];
				$mis=$temp[1];
				$queryinsert=$temp[5];
				$targetinsert=$temp[7];
				$start=$temp[15];
				$end=$temp[16];
				@before=@temp;
			}else {
				if ( ($match+$mis)/($match+$mis+$queryinsert+$targetinsert)<($temp[0]+$temp[1])/($temp[0]+$temp[1]+$temp[5]+$temp[7]) ) {
					$match=$temp[0];
					$mis=$temp[1];
					$queryinsert=$temp[5];
					$targetinsert=$temp[7];
					$start=$temp[15];
					$end=$temp[16];
					@before=@temp;
				}
			}
		}
	}
	print OUT join("\t",@before)."\n";
	close OUT;
	return;
}

sub variant_blat{
	my $blatfile=shift;
	my $targetgenome=$genomefile;
	my $querygenome=$genomenew;
	my $outfile="$prefix.var";
	
	my ($id,%ref,%query,%query_rev,);

	open (IN,"<$targetgenome")||die("fail to open $targetgenome.\n");
	while(<IN>){
		chomp;
		if (/^>(\S+)/) {
			$id=$1;
			next;
		}
		$ref{$id}.=$_;
	}
	close IN;
	foreach $id (keys %ref) {
		$ref_len{$id}=length($ref{$id});
	}

	open (IN,"<$querygenome")||die("fail to open $querygenome.\n");
	while(<IN>){
		chomp;
		if (/^>(\S+)/) {
			$id=$1;
			next;
		}
		$query{$id}.=$_;
	}
	close IN;

	foreach my $id (keys %query) {
		$new_len{$id}=length($query{$id});
		$query_rev{$id}=$query{$id};
		$query_rev{$id}=reverse $query_rev{$id};
		$query_rev{$id}=~tr/AaTtGgCc/TtAaCcGg/;
	}

	my @list=();
	open (IN,"<$blatfile")||die("fail to open $blatfile.\n");
	while (<IN>) {
		chomp;
		#3640	232	0	0	21	240	18	153	+	gi|270610521|gb|GU183365.1|	37320	0	4112	S2846	19749	11158	15183	29	652,720,679,431,276,155,62,32,175,16,8,6,24,18,14,5,8,44,37,106,4,58,157,56,12,51,12,44,10,	0,661,1386,2095,2529,2806,2964,3052,3084,3259,3275,3288,3294,3325,3385,3399,3412,3432,3476,3514,3621,3629,3689,3846,3930,3950,4025,4040,4102,	11158,11810,12533,13212,13645,13921,14081,14144,14182,14360,14380,14388,14397,14421,14439,14456,14464,14510,14555,14599,14705,14714,14772,14952,15039,15065,15116,15129,15173,
		my @temp=split /\t/,$_;
		push @list,\@temp;
	}
	close IN;
	@list=sort {$a->[15] <=> $b->[15]} @list;
	@list=sort {$a->[13] cmp $b->[13]} @list;

	my $tmp=$blatfile."_tmp";
	open (OUT,">$tmp")||die("fail to open $tmp.\n");
	foreach (@list) {
		print OUT join("\t",@$_)."\n";
	}
	close OUT;
	@list=();

	my $poschange="$prefix.pos_change";
	open (POS,">$poschange")||die("fail to open $poschange.\n");
	print POS "#Chr\tOldStart\tOldEnd\tNewStart\tNewEnd\n";

	open (OUT,">$outfile")||die("fail to open $outfile.\n");

	my $query_start=undef;#query end of last blat record 
	my $target_start=undef;#target end of last blat record
	open (IN,"<$tmp")||die("fail to open $tmp.\n");
	while(<IN>){
		chomp;
		#3640	232	0	0	21	240	18	153	+	gi|270610521|gb|GU183365.1|	37320	0	4112	S2846	19749	11158	15183	29	652,720,679,431,276,155,62,32,175,16,8,6,24,18,14,5,8,44,37,106,4,58,157,56,12,51,12,44,10,	0,661,1386,2095,2529,2806,2964,3052,3084,3259,3275,3288,3294,3325,3385,3399,3412,3432,3476,3514,3621,3629,3689,3846,3930,3950,4025,4040,4102,	11158,11810,12533,13212,13645,13921,14081,14144,14182,14360,14380,14388,14397,14421,14439,14456,14464,14510,14555,14599,14705,14714,14772,14952,15039,15065,15116,15129,15173,
		my @temp=split /\t/,$_;
		my @size=split /\,/,$temp[18];
		my @qstart=split /\,/,$temp[19];
		my @tstart=split /\,/,$temp[20];
		my $target=$temp[13];
		my $query=$temp[9];
		if (defined $query_start) {
			next if ($qstart[0]-$query_start<=0);
			if ($temp[8] eq "+") {			
				my $ref_pos1=$target_start;
				my $ref_pos2=$tstart[0];
				my $query_pos1=$query_start;
				my $query_pos2=$qstart[0];
				my $ref_base=substr($ref{$target},$ref_pos1-1,$ref_pos2-$ref_pos1+1);
				my $query_base=substr($query{$query},$query_pos1-1,$query_pos2-$query_pos1+1);
				$ref_base=uc($ref_base);
				$query_base=uc($query_base);
				if ($ref_base ne $query_base) {
					printf OUT "$target\t$ref_pos1\t$query_pos1\t$ref_base\t$query_base\t\.\t\.\tDP=1\n";
					print POS "$temp[13]\t$ref_pos1\t$ref_pos2\t$query_pos1\t$query_pos2\t$temp[9]\t$temp[8]\n";
					$change{$temp[13]}{"$ref_pos1\t$ref_pos2"}="$query_pos1\t$query_pos2\t$temp[9]\t$temp[8]";
				}
			}else {
				my $ref_pos1=$target_start;
				my $ref_pos2=$tstart[0];
				my $query_pos1=$query_start;
				my $query_pos2=$qstart[0];
				my $ref_base=substr($ref{$target},$ref_pos1-1,$ref_pos2-$ref_pos1+1);
				my $query_base=substr($query_rev{$query},$query_pos1-1,$query_pos2-$query_pos1+1);
				$ref_base=uc($ref_base);
				$query_base=uc($query_base);
				if ($ref_base ne $query_base) {
					printf OUT "$target\t$ref_pos1\t$query_pos1\t$ref_base\t$query_base\t\.\t\.\tDP=1\n";
					print POS "$temp[13]\t$ref_pos1\t$ref_pos2\t$query_pos1\t$query_pos2\t$temp[9]\t$temp[8]\n";
					$change{$temp[13]}{"$ref_pos1\t$ref_pos2"}="$query_pos1\t$query_pos2\t$temp[9]\t$temp[8]";
				}
			}
			$target_start=$tstart[$#size]+$size[$#size];
			$query_start=$qstart[$#size]+$size[$#size];
		}
		for (my $j=0;$j<@size;$j++) {
			print POS "$temp[13]\t".($tstart[$j]+1)."\t".($tstart[$j]+$size[$j])."\t".($qstart[$j]+1)."\t".($qstart[$j]+$size[$j])."\t$temp[9]\t$temp[8]\n";
			$change{$temp[13]}{($tstart[$j]+1)."\t".($tstart[$j]+$size[$j])}=($qstart[$j]+1)."\t".($qstart[$j]+$size[$j])."\t$temp[9]\t$temp[8]";
			my $old_start;
			my $new_start;
			for (my $k=1;$k<=$size[$j];$k++) {
				if ($temp[8] eq "+") {
					my $ref_pos=$tstart[$j]+$k;
					$old_start=$ref_pos if (!defined $old_start);
					my $query_pos=$qstart[$j]+$k;
					$new_start=$query_pos if (!defined $new_start);
					my $ref_base=substr($ref{$target},$ref_pos-1,1);
					my $query_base=substr($query{$query},$query_pos-1,1);
					$ref_base=uc($ref_base);
					$query_base=uc($query_base);
					if ($ref_base ne $query_base) {
						printf OUT "$target\t$ref_pos\t$query_pos\t$ref_base\t$query_base\t\.\t\.\tDP=1\n";
					}
				}else {
					my $ref_pos=$tstart[$j]+$k;
					my $query_pos=$qstart[$j]+$k;
					my $ref_base=substr($ref{$target},$ref_pos-1,1);
					my $query_base=substr($query_rev{$query},$query_pos-1,1);
					$ref_base=uc($ref_base);
					$query_base=uc($query_base);
					if ($ref_base ne $query_base) {
						printf OUT "$target\t$ref_pos\t$query_pos\t$ref_base\t$query_base\t\.\t\.\tDP=1\n";
						$change{$temp[13]}{$ref_pos}="$ref_base\t$query_base\t$temp[9]";
					}
				}
			}
		}
		for (my $j=0;$j<@size-1;$j++) {
			if ($temp[8] eq "+") {
				my $ref_pos1=$tstart[$j]+$size[$j];
				my $ref_pos2=$tstart[$j+1];
				my $query_pos1=$qstart[$j]+$size[$j];
				my $query_pos2=$qstart[$j+1];
				my $ref_base=substr($ref{$target},$ref_pos1-1,$ref_pos2-$ref_pos1+1);
				my $query_base=substr($query{$query},$query_pos1-1,$query_pos2-$query_pos1+1);
				$ref_base=uc($ref_base);
				$query_base=uc($query_base);
				if ($ref_base ne $query_base) {
					printf OUT "$target\t$ref_pos1\t$query_pos1\t$ref_base\t$query_base\t\.\t\.\tDP=1\n";
					print POS "$temp[13]\t$ref_pos1\t$ref_pos2\t$query_pos1\t$query_pos2\t$temp[9]\t$temp[8]\n";
					$change{$temp[13]}{"$ref_pos1\t$ref_pos2"}="$query_pos1\t$query_pos2\t$temp[9]\t$temp[8]";
				}
			}else {
				my $ref_pos1=$tstart[$j]+$size[$j];
				my $ref_pos2=$tstart[$j+1];
				my $query_pos1=$qstart[$j]+$size[$j];
				my $query_pos2=$qstart[$j+1];
				my $ref_base=substr($ref{$target},$ref_pos1-1,$ref_pos2-$ref_pos1+1);
				my $query_base=substr($query_rev{$query},$query_pos1-1,$query_pos2-$query_pos1+1);
				$ref_base=uc($ref_base);
				$query_base=uc($query_base);
				if ($ref_base ne $query_base) {
					printf OUT "$target\t$ref_pos1\t$query_pos1\t$ref_base\t$query_base\t\.\t\.\tDP=1\n";
					print POS "$temp[13]\t$ref_pos1\t$ref_pos2\t$query_pos1\t$query_pos2\t$temp[9]\t$temp[8]\n";
					$change{$temp[13]}{"$ref_pos1\t$ref_pos2"}="$query_pos1\t$query_pos2\t$temp[9]\t$temp[8]";
				}
			}
		}
	}
	close IN;
	close OUT;
	close POS;
	return;
}

sub annotation_transfer	{
	my $ref_file=$genomefile;
	my $query_file=$genomenew;

	%ref=();
	%seq=();
	open (IN,"<$ref_file")||die("fail to open $ref_file.\n");
	my $id=undef;
	while (<IN>) {
		chomp;
		if (/^>(\S+)/) {
			$id=$1;
			next;
		}
		$ref{$id}.=$_;
	}
	close IN;

	open (IN,"<$query_file")||die("fail to open $query_file.\n");
	$id=undef;
	while (<IN>) {
		chomp;
		if (/^>(\S+)/) {
			$id=$1;
			next;
		}
		$seq{$id}.=$_;
	}
	close IN;

	foreach my $id (keys %change) {
		foreach my $pos (keys %{$change{$id}}) {
			my ($ref_pos1,$ref_pos2)=split /\t/,$pos;
			my ($query_pos1,$query_pos2,$chr,$strand)=split /\t/,$change{$id}{$pos};
			my $num=$ref_pos2-$ref_pos1+1;
			$num=$query_pos2-$query_pos1+1 if ($num>$query_pos2-$query_pos1+1);
			for (my $j=0;$j<$num;$j++) {
				if ($strand eq "+") {
					$pos{$id}{$ref_pos1+$j}=($query_pos1+$j)."\t$chr";
				}else {
					$pos{$id}{$ref_pos1+$j}=(abs(length($seq{$chr})-($query_pos1+$j))+1)."\t$chr";
				}
			}
		}
	}

	my $annout="$prefix.anno_update";
	if ($annfile=~/\.tbl$/) {
		&tblupdate2($annfile,"$annout");
	}elsif ($annfile=~/\.gtf$/ or $annfile=~/\.gff$/ or $annfile=~/\.gff2$/ or $annfile=~/\.gff3$/) {
		&gffupdate2($annfile,"$annout");
	}elsif ($annfile=~/\.bed$/) {
		&bedupdate2($annfile,"$annout");
	}
	return;
}

sub tblupdate2 {
	my $infile=shift;
	my $outfile=shift;

	my $outfile1=$outfile;
	my $outfile2=$outfile.".check";
	open (OUT1,">$outfile1")||die("fail to open $outfile1.\n");
	open (OUT2,">$outfile2")||die("fail to open $outfile2.\n");
	my $id=undef;
	my $query=undef;
	my $type=undef;
	my $typenew=undef;
	my $strand=undef;
	my $length=undef;
	my @str=();
	my %type1=();
	my %type2=();
	my %type3=();
	my $block=0;
	my $flag=0;
	my $stop_codon=0;
	my $query_pos1;
	my $query_pos2;
	my $query_chr1;
	my $query_chr2;
	open (IN,"<$infile")||die("fail to open $infile.\n");
	while (<IN>) {
		chomp;
		if (/^>\S+\s+(\S+)$/) {
			#>Feature gi|444893469|emb|AL123456.3|
			$id=$1;
		}elsif (/^<?(\d+)\t(\d+)\t(\S+)$/) {
			#122026460	125184587	centromere
			#11874	14409	gene
			$typenew=$3;
			$stop_codon=0;
			if ($1<$2) {
				my $seq=substr($ref{$id},$1-1,$2-$1+1);
				for (my $j=0;$j<length($seq);$j+=3) {
					my $codon=substr($seq,$j,3);
					$stop_codon++ if ((&codon2amino($codon)) eq "*");
				}
			}else {
				my $seq=substr($ref{$id},$2-1,$1-$2+1);
				$seq=reverse $seq;
				$seq=~tr/AaTtGgCc/TtAaCcGg/;
				for (my $j=0;$j<length($seq);$j+=3) {
					my $codon=substr($seq,$j,3);
					$stop_codon++ if ((&codon2amino($codon)) eq "*");
				}
			}
			$length=abs($2-$1)+1;
			$block=0;
			$flag=0;
			$type1{$typenew}++;
			$query_pos1="?";$query_chr1="?";
			$query_pos2="?";$query_chr2="?";
			($query_pos1,$query_chr1)=split /\t/,$pos{$id}{$1};
			($query_pos2,$query_chr2)=split /\t/,$pos{$id}{$2};
			if ($query_chr1 ne "?" and $query_chr1 eq $query_chr2) {
				if (abs($query_pos2-$query_pos1)+1>2*$length) {
					$query_pos2="?";
				}
			}elsif ($query_chr1 ne "?" and $query_chr2 ne "?" and $query_chr1 ne $query_chr2) {
				$query_pos2="?";
			}
			if ($query_pos1 ne "?" and $query_pos2 ne "?") {
				$flag=1;
				$block++;
				if ($query_pos1<$query_pos2) {
					$strand="+";
				}else {
					$strand="-";
				}
				if (!defined $query or $query ne $query_chr1) {
					$query=$query_chr1;
					print OUT1 ">Feature $query\n";
					print OUT2 ">Feature $query\n";
				}
			}elsif ($query_pos1 ne "?") {
				$block++;
				my $query_pos="?";my $query_chr="?";
				if ($1<$2) {
					($query_pos,$query_chr)=split /\t/,$pos{$id}{$1+1};
					if ($query_pos ne "?") {
						$flag=1;
						if ($query_pos1<$query_pos) {
							$strand="+";
						}else {
							$strand="-";
						}
					}
				}else {
					($query_pos,$query_chr)=split /\t/,$pos{$id}{$1-1};
					if ($query_pos ne "?") {
						$flag=1;
						if ($query_pos1<$query_pos) {
							$strand="+";
						}else {
							$strand="-";
						}
					}
				}
				if (!defined $query or $query ne $query_chr1) {
					$query=$query_chr1;
					print OUT1 ">Feature $query\n";
					print OUT2 ">Feature $query\n";
				}
			}elsif ($query_pos2 ne "?") {
				$block++;
				my $query_pos="?";my $query_chr="?";
				if ($1<$2) {
					($query_pos,$query_chr)=split /\t/,$pos{$id}{$2-1};
					if ($query_pos ne "?") {
						$flag=1;
						if ($query_pos<$query_pos2) {
							$strand="+";
						}else {
							$strand="-";
						}
					}
				}else {
					($query_pos,$query_chr)=split /\t/,$pos{$id}{$2+1};
					if ($query_pos ne "?") {
						$flag=1;
						if ($query_pos<$query_pos2) {
							$strand="+";
						}else {
							$strand="-";
						}
					}
				}
				if (!defined $query or $query ne $query_chr2) {
					$query=$query_chr2;
					print OUT1 ">Feature $query\n";
					print OUT2 ">Feature $query\n";
				}
			}
			if (@str>0) {
				if (@str==$block) {
					my ($translate,$print_str)=split /\:/,(&block_print(\@str,$type,$query));
					if (defined $translate and $translate==0) {
						$type2{$type}++;
						print OUT1 "$print_str";
					}elsif (defined $translate and $translate==1) {
						$type2{$type}++;
						$type3{$type}++;
						print OUT1 "$print_str";
						print OUT2 "$print_str\t\t\tnote\tTerminated by stop codon\n";
					}
				}elsif (@str<$block) {
					my ($translate,$print_str)=split /\:/,(&block_print(\@str,$type,$query));
					if (defined $print_str) {
						$type2{$type}++;
						$type3{$type}++;
						print OUT1 "$print_str";
						print OUT2 "$print_str\t\t\tnote\tPartial tranferred\n";
					}
				}
			}
			@str=();
			$type=$typenew;
			if ($flag==1) {
				push (@str,"$query_pos1\t$query_pos2\t$strand\t$length\t$stop_codon");
			}
		}elsif (/^\t\t\t\S+\t.+$/) {
			#			product	DEAD/H (Asp-Glu-Ala-Asp/His) box helicase 11 like 1
			#			note	Linear centromere model derived predominantly from reads generated in PMID: 17803354. This region does not represent an actual centromere sequence, as long-range ordering of repeats and unmapped WGS contigs is not provided by the model. For details of model production, see http://arxiv.org/abs/1307.0035.
			#			gene	DDX11L1
			if (@str>0) {
				if (@str==$block) {
					my ($translate,$print_str)=split /\:/,(&block_print(\@str,$type,$query));
					if (defined $translate and $translate==0) {
						$type2{$type}++;
						print OUT1 "$print_str";
					}elsif (defined $translate and $translate==1) {
						$type2{$type}++;
						$type3{$type}++;
						print OUT1 "$print_str";
						print OUT2 "$print_str\t\t\tnote\tTerminated by stop codon\n";
					}
				}elsif (@str<$block) {
					my ($translate,$print_str)=split /\:/,(&block_print(\@str,$type,$query));
					if (defined $print_str) {
						$type2{$type}++;
						$type3{$type}++;
						print OUT1 "$print_str";
						print OUT2 "$print_str\t\t\tnote\tPartial tranferred\n";
					}
				}
			}
			@str=();
			if ($flag==1) {
				print OUT1 "$_\n";
			}
		}elsif (/^<?(\d+)\t(\d+)$/) {
			$length=abs($2-$1)+1;
			$query_pos1="?";$query_chr1="?";
			$query_pos2="?";$query_chr2="?";
			($query_pos1,$query_chr1)=split /\t/,$pos{$id}{$1};
			($query_pos2,$query_chr2)=split /\t/,$pos{$id}{$2};
			if ($query_chr1 eq $query_chr2) {
				if (abs($query_pos2-$query_pos1)>2*$length) {
					$query_pos2="?";
				}
			}else {
				$query_pos2="?";
			}
			if ($query_pos1 ne "?" and $query_pos2 ne "?") {
				$flag=1;
				$block++;
				if ($query_pos1<$query_pos2) {
					$strand="+";
				}else {
					$strand="-";
				}
				if (!defined $query or $query ne $query_chr1) {
					$query=$query_chr1;
					print OUT1 ">Feature $query\n";
					print OUT2 ">Feature $query\n";
				}
			}elsif ($query_pos1 ne "?") {
				$block++;
				my $query_pos="?";my $query_chr="?";
				if ($1<$2) {
					($query_pos,$query_chr)=split /\t/,$pos{$id}{$1+1};
					if ($query_pos ne "?") {
						$flag=1;
						if ($query_pos1<$query_pos) {
							$strand="+";
						}else {
							$strand="-";
						}
					}
				}else {
					($query_pos,$query_chr)=split /\t/,$pos{$id}{$1-1};
					if ($query_pos ne "?") {
						$flag=1;
						if ($query_pos1<$query_pos) {
							$strand="+";
						}else {
							$strand="-";
						}
					}
				}
				if (!defined $query or $query ne $query_chr1) {
					$query=$query_chr1;
					print OUT1 ">Feature $query\n";
					print OUT2 ">Feature $query\n";
				}
			}elsif ($query_pos2 ne "?") {
				$block++;
				my $query_pos="?";my $query_chr="?";
				if ($1<$2) {
					($query_pos,$query_chr)=split /\t/,$pos{$id}{$2-1};
					if ($query_pos ne "?") {
						$flag=1;
						if ($query_pos<$query_pos2) {
							$strand="+";
						}else {
							$strand="-";
						}
					}
				}else {
					($query_pos,$query_chr)=split /\t/,$pos{$id}{$2+1};
					if ($query_pos ne "?") {
						$flag=1;
						if ($query_pos<$query_pos2) {
							$strand="+";
						}else {
							$strand="-";
						}
					}
				}
				if (!defined $query or $query ne $query_chr2) {
					$query=$query_chr2;
					print OUT1 ">Feature $query\n";
					print OUT2 ">Feature $query\n";
				}
			}
			if (exists $pos{$id}{$1} or exists $pos{$id}{$2}) {
				push (@str,"$query_pos1\t$query_pos2\t$strand\t$length\t0\n");
			}
		}else {
			if ($flag==1) {
				print OUT1 "$_\n";
			}
		}
	}
	close IN;
	if (@str>0) {
		if (@str==$block) {
			my ($translate,$print_str)=split /\:/,(&block_print(\@str,$type,$query));
			if (defined $translate and $translate==0) {
				$type2{$type}++;
				print OUT1 "$print_str";
			}elsif (defined $translate and $translate==1) {
				$type2{$type}++;
				$type3{$type}++;
				print OUT1 "$print_str";
				print OUT2 "$print_str\t\t\tnote\tTerminated by stop codon\n";
			}
		}elsif (@str<$block) {
			my ($translate,$print_str)=split /\:/,(&block_print(\@str,$type,$query));
			if (defined $print_str) {
				$type2{$type}++;
				$type3{$type}++;
				print OUT1 "$print_str";
				print OUT2 "$print_str\t\t\tnote\tPartial tranferred\n";
			}
		}
	}
	@str=();
	close OUT1;
	close OUT2;
	print "#Raw element:\n";
	foreach my $type (sort (keys %type1)) {
		print "$type\t$type1{$type}\n";
	}
	print "#Complete element:\n";
	foreach my $type (sort (keys %type2)) {
		print "$type\t$type2{$type}\n";
	}
	print "#Check element:\n";
	foreach my $type (sort (keys %type3)) {
		print "$type\t$type3{$type}\n";
	}
	return;
}

sub block_print {
	my $str=shift;
	my @str=@{$str};
	my $type=shift;
	my $query=shift;

	my $print_str=undef;
	my $extend_type=0;
	my $extend=300;
	my $translate=0;
	for (my $i=0;$i<@str;$i++) {
		my ($pos1,$pos2,$strand,$length,$stop_codon)=split /\t/,$str[$i];
		my $flag=0;
		if ($pos1 ne "" and $pos2 ne "") {
			$flag=0;
		}elsif ($pos1 ne "") {
			$flag=1;
			if ($strand eq "+") {
				if ($pos1+$length-1<length($seq{$query})) {
					$pos2=$pos1+$length-1;
				}else {
					$pos2=length($seq{$query});
				}
			}else {
				if ($pos1-$length+1>0) {
					$pos2=$pos1-$length+1;
				}else {
					$pos2=1;
				}
			}
		}elsif ($pos2 ne "") {
			$flag=2;
			if ($strand eq "+") {
				if ($pos2-$length+1>0) {
					$pos1=$pos2-$length+1;
				}else {
					$pos1=1;
				}
			}else {
				if ($pos2+$length-1<length($seq{$query})) {
					$pos1=$pos2+$length-1;
				}else {
					$pos1=length($seq{$query});
				}
			}
		}
		if ($i==0) {
			if ($pos1<$pos2) {
				if ($type eq "CDS") {
					#(int(($length/2)/3))*3
					if ($flag==1) {
						my $seq=substr($seq{$query},$pos1-1,$pos2-$pos1+1);
						my %codon=();
						for (my $j=0;$j<3;$j++) {
							$codon{$j}=0;
							for (my $k=$j;$k<length($seq)-3;$k+=3) {
								my $codon=substr($seq,$k,3);
								$codon{$j}++ if ((&codon2amino($codon)) eq "*");
							}
						}
						my @codon=sort {$codon{$a} <=> $codon{$b}} keys %codon;
						$pos1=$pos1+$codon[0];
						$pos2=$pos2+$codon[0];
						if ($codon{$codon[0]}>$stop_codon) {
							$translate=1;
						}
						for (my $j=$pos1+(int(($length/2)/3))*3-1;$j>$pos1-$extend and $j>=0;$j-=3) {
							my $codon=substr($seq{$query},$j,3);
							my $amino=&codon2amino($codon);
							if ($amino eq "M") {
								if ($j+1-$pos1>=$length*0.2) {
								}else {
									$pos1=$j+1;
									last;
								}
							}elsif ($amino eq "*") {
								if ($j+4-$pos1>=$length*0.2) {
									last;
								}else {
									$pos1=$j+4;
									last;
								}
							}
						}
						if (@str==1) {
							for (my $j=$pos1+(int(($length/2)/3))*3-1;$j<$pos2+$extend;$j+=3) {
								my $codon=substr($seq{$query},$j,3);
								my $amino=&codon2amino($codon);
								if ($amino eq "*") {
									if ($pos2-($j+3)>=$length*0.2) {
										last;
									}else {
										$pos2=$j+3;
										last;
									}
								}
							}
						}
					}elsif ($flag==2) {
						my $seq=substr($seq{$query},$pos1-1,$pos2-$pos1+1);
						my %codon=();
						for (my $j=0;$j<3;$j++) {
							$codon{$j}=0;
							for (my $k=length($seq)-6+$j;$k>=0;$k-=3) {
								my $codon=substr($seq,$k,3);
								$codon{$j}++ if ((&codon2amino($codon)) eq "*");
							}
						}
						my @codon=sort {$codon{$a} <=> $codon{$b}} keys %codon;
						$pos1=$pos1+$codon[0];
						$pos2=$pos2+$codon[0];
						if ($codon{$codon[0]}>$stop_codon) {
							$translate=1;
						}
						for (my $j=$pos2-(int(($length/2)/3))*3-3;$j>$pos1-$extend and $j>=0;$j-=3) {
							my $codon=substr($seq{$query},$j,3);
							my $amino=&codon2amino($codon);
							if ($amino eq "M") {
								if ($j+1-$pos1>=$length*0.2) {
								}else {
									$pos1=$j+1;
									last;
								}
							}elsif ($amino eq "*") {
								if ($j+4-$pos1>=$length*0.2) {
									last;
								}else {
									$pos1=$j+4;
									last;
								}
							}
						}
						if (@str==1) {
							for (my $j=$pos2-(int(($length/2)/3))*3-3;$j<$pos2+$extend;$j+=3) {
								my $codon=substr($seq{$query},$j,3);
								my $amino=&codon2amino($codon);
								if ($amino eq "*") {
									if ($pos2-($j+3)>=$length*0.2) {
										last;
									}else {
										$pos2=$j+3;
										last;
									}
								}
							}
						}
					}
				}
			}else {
				if ($type eq "CDS") {
					if ($flag==1) {
						my $seq=substr($seq{$query},$pos2-1,$pos1-$pos2+1);
						$seq=reverse $seq;
						$seq=~tr/AaTtGgCc/TtAaCcGg/;
						my %codon=();
						for (my $j=0;$j<3;$j++) {
							$codon{$j}=0;
							for (my $k=$j;$k<length($seq)-3;$k+=3) {
								my $codon=substr($seq,$k,3);
								$codon{$j}++ if ((&codon2amino($codon)) eq "*");
							}
						}
						my @codon=sort {$codon{$a} <=> $codon{$b}} keys %codon;
						$pos1=$pos1-$codon[0];
						$pos2=$pos2-$codon[0];
						if ($codon{$codon[0]}>$stop_codon) {
							$translate=1;
						}
						for (my $j=$pos1-(int(($length/2)/3))*3-3;$j<$pos1+$extend;$j+=3) {
							my $codon=substr($seq{$query},$j,3);
							$codon=reverse $codon;
							$codon=~tr/AaTtGgCc/TtAaCcGg/;
							my $amino=&codon2amino($codon);
							if ($amino eq "M") {
								if ($pos1-($j+3)>=$length*0.2) {
								}else {
									$pos1=$j+3;
									last;
								}
							}elsif ($amino eq "*") {
								if ($pos1-$j>=$length*0.2) {
									last;
								}else {
									$pos1=$j;
									last;
								}
							}
						}
						if (@str==1) {
							for (my $j=$pos1-(int(($length/2)/3))*3-3;$j>$pos2-$extend and $j>=0;$j-=3) {
								my $codon=substr($seq{$query},$j,3);
								$codon=reverse $codon;
								$codon=~tr/AaTtGgCc/TtAaCcGg/;
								my $amino=&codon2amino($codon);
								if ($amino eq "*") {
									if ($j+1-$pos2>=$length*0.2) {
										last;
									}else {
										$pos2=$j+1;
										last;
									}
								}
							}
						}
					}elsif ($flag==2) {
						my $seq=substr($seq{$query},$pos2-1,$pos1-$pos2+1);
						$seq=reverse $seq;
						$seq=~tr/AaTtGgCc/TtAaCcGg/;
						my %codon=();
						for (my $j=0;$j<3;$j++) {
							$codon{$j}=0;
							for (my $k=length($seq)-6+$j;$k>=3;$k-=3) {
								my $codon=substr($seq,$k,3);
								$codon{$j}++ if ((&codon2amino($codon)) eq "*");
							}
						}
						my @codon=sort {$codon{$a} <=> $codon{$b}} keys %codon;
						$pos1=$pos1-$codon[0];
						$pos2=$pos2-$codon[0];
						if ($codon{$codon[0]}>$stop_codon) {
							$translate=1;
						}
						for (my $j=$pos2+(int(($length/2)/3))*3-1;$j<$pos1+$extend;$j+=3) {
							my $codon=substr($seq{$query},$j,3);
							$codon=reverse $codon;
							$codon=~tr/AaTtGgCc/TtAaCcGg/;
							my $amino=&codon2amino($codon);
							if ($amino eq "M") {
								if ($pos1-($j+3)>=$length*0.2) {
								}else {
									$pos1=$j+3;
									last;
								}
							}elsif ($amino eq "*") {
								if ($pos1-$j>=$length*0.2) {
									last;
								}else {
									$pos1=$j;
									last;
								}
							}
						}
						if (@str==1) {
							for (my $j=$pos2+(int(($length/2)/3))*3-1;$j>$pos2-$extend and $j>=0;$j-=3) {
								my $codon=substr($seq{$query},$j,3);
								$codon=reverse $codon;
								$codon=~tr/AaTtGgCc/TtAaCcGg/;
								my $amino=&codon2amino($codon);
								if ($amino eq "*") {
									if ($j+1-$pos2>=$length*0.2) {
										last;
									}else {
										$pos2=$j+1;
										last;
									}
								}
							}
						}
					}
				}
			}
			$print_str.="$pos1\t$pos2\t$type\n";
		}elsif ($i>0 and $i<@str-1) {
			$print_str.="$pos1\t$pos2\n";
		}else {
			if ($pos1<$pos2) {
				if ($type eq "CDS") {
					if ($flag==1) {
						for (my $j=$pos1+(int(($length/2)/3))*3-1;$j<$pos2+$extend;$j+=3) {
							my $codon=substr($seq{$query},$j,3);
							my $amino=&codon2amino($codon);
							if ($amino eq "*") {
								if ($pos2-($j+3)>=$length*0.2) {
									last;
								}else {
									$pos2=$j+3;
									last;
								}
							}
						}
					}elsif ($flag==2) {
						for (my $j=$pos2-(int(($length/2)/3))*3-3;$j<$pos2+$extend;$j+=3) {
							my $codon=substr($seq{$query},$j,3);
							my $amino=&codon2amino($codon);
							if ($amino eq "*") {
								if ($pos2-($j+3)>=$length*0.2) {
									last;
								}else {
									$pos2=$j+3;
									last;
								}
							}
						}
					}
				}
			}else {
				if ($type eq "CDS") {
					if ($flag==1) {
						for (my $j=$pos1-(int(($length/2)/3))*3-3;$j>$pos2-$extend and $j>=0;$j-=3) {
							my $codon=substr($seq{$query},$j,3);
							$codon=reverse $codon;
							$codon=~tr/AaTtGgCc/TtAaCcGg/;
							my $amino=&codon2amino($codon);
							if ($amino eq "*") {
								if ($j+1-$pos2>=$length*0.2) {
									last;
								}else {
									$pos2=$j+1;
									last;
								}
							}
						}
					}elsif ($flag==2) {
						for (my $j=$pos2+(int(($length/2)/3))*3-1;$j>$pos2-$extend and $j>=0;$j-=3) {
							my $codon=substr($seq{$query},$j,3);
							$codon=reverse $codon;
							$codon=~tr/AaTtGgCc/TtAaCcGg/;
							my $amino=&codon2amino($codon);
							if ($amino eq "*") {
								if ($j+1-$pos2>=$length*0.2) {
									last;
								}else {
									$pos2=$j+1;
									last;
								}
							}
						}
					}
				}
			}
			$print_str.="$pos1\t$pos2\n";
		}
	}
	return "$translate\:$print_str";
}

sub gffupdate2 {
	my $infile=shift;
	my $outfile=shift;
	
	open (OUT,">$outfile")||die("fail to open $outfile.\n");
	open (IN,"<$infile")||die("fail to open $infile.\n");
	while (<IN>) {
		chomp;
		#chr1	Cufflinks	exon	4886744	4886831	.	+	.	gene_id "XLOC_000005"; transcript_id "TCONS_00000014"; exon_number "4"; oId "CUFF.13.9"; tss_id "TSS8";
		if ($_=~/^\#/) {
			next;
		}
		my @list=split/\t/,$_;
		my $query_pos1="?";my $query_chr1="?";
		my $query_pos2="?";my $query_chr2="?";
		($query_pos1,$query_chr1)=split /\t/,$pos{$list[0]}{$list[3]};
		($query_pos2,$query_chr2)=split /\t/,$pos{$list[0]}{$list[4]};
		if ($query_chr1 eq $query_chr2) {
			$list[0]=$query_chr1;
			if ($query_pos1<=$query_pos2) {
				$list[3]=$query_pos1;
				$list[4]=$query_pos2;
				$list[6]="+";
			}elsif ($query_pos1>$query_pos2) {
				$list[3]=$query_pos2;
				$list[4]=$query_pos1;
				$list[6]="-";
			}
			my $list=join "\t",@list;
			print OUT "$list\n";
		}
	}
	close(IN);
	close OUT;
	return;
}

sub bedupdate2 {
	my $infile=shift;
	my $outfile=shift;
	
	open (OUT,">$outfile")||die("fail to open $outfile.\n");
	open (IN,"<$infile")||die("fail to open $infile.\n");
	while (<IN>) {
		chomp;
		#chr22 1000 5000 cloneA 960 + 1000 5000 0 2 567,488, 0,3512
		if ($_=~/^\#/ or $_=~/^track/) {
			next;
		}
		my @list=split/\s+/,$_;
		my $query_pos1="?";
		my $query_chr1="?";
		($query_pos1,$query_chr1)=split /\t/,$pos{$list[0]}{$list[1]+1};
		my $query_pos2="?";
		my $query_chr2="?";
		($query_pos2,$query_chr2)=split /\t/,$pos{$list[0]}{$list[2]};
		my $len;
		my $sta;
		if ($query_chr1 ne "?" and $query_chr1 eq $query_chr2) {
			if ($query_pos1<=$query_pos2) {
				if (defined $list[11]) {
					$len="";
					$sta="";
					my @len=split /\,/,$list[10];
					my @start=split /\,/,$list[11];
					for (my $i=0;$i<@start;$i++) {
						my $len_new="?";
						my $sta_new="?";
						$sta_new=$pos{$list[0]}{$list[1]+$start[$i]+1}-$pos{$list[0]}{$list[1]+1} if (exists $pos{$list[0]}{$list[1]+1} and exists $pos{$list[0]}{$list[1]+$start[$i]+1});
						$len_new=$pos{$list[0]}{$list[1]+$start[$i]+$len[$i]}-$pos{$list[0]}{$list[1]+$start[$i]+1}+1 if (exists $pos{$list[0]}{$list[1]+$start[$i]+$len[$i]} and exists $pos{$list[0]}{$list[1]+$start[$i]+1});
						$sta.=$sta_new.",";
						$len.=$len_new.",";
					}
				}
				$list[0]=$query_chr1;
				$list[1]=$query_pos1-1;
				$list[2]=$query_pos2;
				$list[5]="+" if (defined $list[5]);
				$list[6]=$query_pos1-1 if (defined $list[6]);
				$list[7]=$query_pos2 if (defined $list[7]);
				@list=(@list[0..7],@list[8..9],$len,$sta);
			}elsif ($query_pos1>$query_pos2) {
				if (defined $list[11]) {
					$len="";
					$sta="";
					my @len=split /\,/,$list[10];
					my @start=split /\,/,$list[11];
					for (my $i=0;$i<@start;$i++) {
						my $len_new="?";
						my $sta_new="?";
						$len_new=abs($pos{$list[0]}{$list[1]+$start[$i]+$len[$i]}-$pos{$list[0]}{$list[1]+$start[$i]+1})+1 if (exists $pos{$list[0]}{$list[1]+$start[$i]+$len[$i]} and exists $pos{$list[0]}{$list[1]+$start[$i]+1});
						$sta_new=abs($pos{$list[0]}{$list[1]+$start[$i]+1}-$len_new+1-$pos{$list[0]}{$list[2]}) if (exists $pos{$list[0]}{$list[2]} and exists $pos{$list[0]}{$list[1]+$start[$i]+1});
						$sta=$sta_new.",".$sta;
						$len=$len_new.",".$len;
					}
				}
				$list[1]=$query_pos2-1;
				$list[2]=$query_pos1;
				$list[5]="-" if (defined $list[5]);
				$list[6]=$query_pos2-1 if (defined $list[6]);
				$list[7]=$query_pos1 if (defined $list[7]);
				@list=(@list[0..7],@list[8..9],$len,$sta);
			}
			my $list=join "\t",@list;
			print OUT "$list\n";
		}
	}
	close(IN);
	close OUT;
	return;
}

sub codon2amino {
	my $codon=$_[0];
	$codon = uc($codon);
	my(%genetic_code)=(
	'TCA' => 'S', # Serine ££Àø∞±À· 
	'TCC' => 'S', # Serine ££Àø∞±À· 
	'TCG' => 'S', # Serine ££Àø∞±À· 
	'TCT' => 'S', # Serine ££Àø∞±À· 
	'TTC' => 'F', # Phenylalanine ££±Ω±˚∞±À· 
	'TTT' => 'F', # Phenylalanine ££±Ω±˚∞±À· 
	'TTA' => 'L', # Leucine ££¡¡∞±À· 
	'TTG' => 'L', # Leucine ££¡¡∞±À· 
	'TAC' => 'Y', # Tyrosine ££¿“∞±À· 
	'TAT' => 'Y', # Tyrosine ££¿“∞±À· 
	'TAA' => '*', # Stop ££Õ£÷π 
	'TAG' => '*', # Stop ££Õ£÷π 
	'TGC' => 'C', # Cysteine ££∞ÎÎ◊∞±À· 
	'TGT' => 'C', # Cysteine ££∞ÎÎ◊∞±À· 
	'TGA' => '*', # Stop ££Õ£÷π 
	'TGG' => 'W', # Tryptophan ££…´∞±À· 
	'CTA' => 'L', # Leucine ££¡¡∞±À· 
	'CTC' => 'L', # Leucine ££¡¡∞±À· 
	'CTG' => 'L', # Leucine ££¡¡∞±À· 
	'CTT' => 'L', # Leucine ££¡¡∞±À· 
	'CCA' => 'P', # Proline ££∏¨∞±À· 
	'CAT' => 'H', # Histidine ££◊È∞±À· 
	'CAA' => 'Q', # Glutamine ££π»∞±ı£∞∑ 
	'CAG' => 'Q', # Glutamine ££π»∞±ı£∞∑ 
	'CGA' => 'R', # Arginine ££æ´∞±À· 
	'CGC' => 'R', # Arginine ££æ´∞±À· 
	'CGG' => 'R', # Arginine ££æ´∞±À· 
	'CGT' => 'R', # Arginine ££æ´∞±À· 
	'ATA' => 'I', # Isoleucine ££“Ï¡¡∞±À· 
	'ATC' => 'I', # Isoleucine ££“Ï¡¡∞±À· 
	'ATT' => 'I', # Isoleucine ££“Ï¡¡∞±À· 
	'ATG' => 'M', # Methionine ££µ∞∞±À· 
	'ACA' => 'T', # Threonine ££À’∞±À· 
	'ACC' => 'T', # Threonine ££À’∞±À· 
	'ACG' => 'T', # Threonine ££À’∞±À· 
	'ACT' => 'T', # Threonine ££À’∞±À· 
	'AAC' => 'N', # Asparagine ££ÃÏ∂¨ı£∞∑ 
	'AAT' => 'N', # Asparagine ££ÃÏ∂¨ı£∞∑ 
	'AAA' => 'K', # Lysine ££¿µ∞±À· 
	'AAG' => 'K', # Lysine ££¿µ∞±À· 
	'AGC' => 'S', # Serine ££Àø∞±À· 
	'AGT' => 'S', # Serine ££Àø∞±À· 
	'AGA' => 'R', # Arginine ££æ´∞±À· 
	'AGG' => 'R', # Arginine ££æ´∞±À· 
	'CCC' => 'P', # Proline ££∏¨∞±À· 
	'CCG' => 'P', # Proline ££∏¨∞±À· 
	'CCT' => 'P', # Proline ££∏¨∞±À· 
	'CAC' => 'H', # Histidine ££◊È∞±À· 
	'GTA' => 'V', # Valine ££Á”∞±À· 
	'GTC' => 'V', # Valine ££Á”∞±À· 
	'GTG' => 'V', # Valine ££Á”∞±À· 
	'GTT' => 'V', # Valine ££Á”∞±À· 
	'GCA' => 'A', # Alanine ££±˚∞±À· 
	'GCC' => 'A', # Alanine ££±˚∞±À· 
	'GCG' => 'A', # Alanine ££±˚∞±À· 
	'GCT' => 'A', # Alanine ££±˚∞±À· 
	'GAC' => 'D', # Aspartic Acid ££ÃÏ∂¨∞±À· 
	'GAT' => 'D', # Aspartic Acid ££ÃÏ∂¨∞±À· 
	'GAA' => 'E', # Glutamic Acid ££π»∞±À· 
	'GAG' => 'E', # Glutamic Acid ££π»∞±À· 
	'GGA' => 'G', # Glycine ££∏ ∞±À· 
	'GGC' => 'G', # Glycine ££∏ ∞±À· 
	'GGG' => 'G', # Glycine ££∏ ∞±À· 
	'GGT' => 'G' # Glycine ££∏ ∞±À· 
	);
	return $genetic_code{$codon};
}

sub genome_compare_figure {
	my (%qgene,%tgene,);

	# width of figure
	my $width = 3000;
	# height of figure
	my $height = 500;
	my $margin = 100;
	my $svg;
	my $xpels1;
	my $xpels2;
	my $xpels;
	my $query_len=undef;
	my $target_len=undef;
	my $query=undef;
	my $target=undef;
	my $strand=undef;

	#read annotation
	if (exists $opts{a}) {
		if ($annfile=~/\.tbl$/) {
			my $tgene=&tbl2exon($annfile);
			%tgene=%{$tgene};
			my $qgene=&tbl2exon("$prefix.anno_update");
			%qgene=%{$qgene};
		}elsif ($annfile=~/\.gtf$/) {
			my $tgene=&gff2exon($annfile);
			%tgene=%{$tgene};
			my $qgene=&gff2exon("$prefix.anno_update");
			%qgene=%{$qgene};
		}elsif ($annfile=~/\.gff$|\.gff2$|\.gff3$/) {
			my $tgene=&gff2exon($annfile);
			%tgene=%{$tgene};
			my $qgene=&gff2exon("$prefix.anno_update");
			%qgene=%{$qgene};
		}elsif ($annfile=~/\.bed$/) {
			my $tgene=&bed2exon($annfile);
			%tgene=%{$tgene};
			my $qgene=&bed2exon("$prefix.anno_update");
			%qgene=%{$qgene};
		}
	}

	open (IN,"<$prefix.pos_change")||die("fail to open $prefix.pos_change.\n");
	while (<IN>) {
		next if (/^\#/);
		chomp;
		#"$id\t$old_start\t$old_end\t$new_start\t$new_end\tNewGenome\t+\n";
		my @list=split /\t/,$_;
		if (!defined $target or $target ne $list[0]) {
			if (defined $target) {
				print O $svg->xmlify();
				close O;
			}
			$target=$list[0];
			$target_len=$ref_len{$target};
			if ($list[5] ne "NewGenome") {
				$query=$list[5];
				$query_len=$new_len{$query};
			}else {
				$query=$list[0];
				$query_len=$new_len{$query};
			}
			$strand=$list[-1];
			open (O,">$prefix.$target\_$query.svg") || die "$!";
			$svg = SVG->new('width',$width,'height',$height);
			$xpels1=($width-2*$margin)/$target_len;
			$xpels2=($width-2*$margin)/$query_len;
			$svg->text('x',$width/2,'y',20,'font-family','TimesNewRoman','font-size',30,'font-color',"black",'text-anchor','middle')->cdata("$target");
			$svg->rect('x',$margin,'y',$margin-20,'width',$width-2*$margin,'height',20,'fill',"white",'stroke',"black",'stroke-width',1);
			foreach my $pos (keys %tgene) {
				my ($start,$end)=split /\t/,$pos;
				my $color=$tgene{$pos}[0];
				my $name=$tgene{$pos}[1];
				$svg->text('x',$margin+$start*$xpels1+($end-$start+1)/2*$xpels1,'y',$margin-35,'font-family','TimesNewRoman','font-size',10,'font-color',"$color",'writing-mode','tb','text-anchor','end')->cdata("$name");
				$svg->rect('x',$margin+$start*$xpels1,'y',$margin-20,'width',($end-$start+1)*$xpels1,'height',20,'fill',"$color",'stroke',"$color",'stroke-width',1);
			}
			$svg->text('x',$margin,'y',$margin-25,'font-family','TimesNewRoman','font-size',10,'font-color',"black",'text-anchor','end')->cdata("0");
			$svg->line('x1',$margin,'y1',$margin-20,'x2',$margin,'y2',$margin-25,'stroke',"black",'stroke-width',1);
			$svg->text('x',$width-$margin,'y',$margin-25,'font-family','TimesNewRoman','font-size',10,'font-color',"black",'text-anchor','start')->cdata( (sprintf "%.2e",$target_len) );
			$svg->line('x1',$width-$margin,'y1',$margin-20,'x2',$width-$margin,'y2',$margin-25,'stroke',"black",'stroke-width',1);
			for (my $i=1;$i<$target_len/20;$i++) {
				$svg->text('x',$margin+$i*($width-2*$margin)/20,'y',$margin-25,'font-family','TimesNewRoman','font-size',10,'font-color',"black",'text-anchor','middle')->cdata( (sprintf "%.2e",$i*($width-2*$margin)/20/$xpels1) );
				$svg->line('x1',$margin+$i*($width-2*$margin)/20,'y1',$margin-20,'x2',$margin+$i*($width-2*$margin)/20,'y2',$margin-25,'stroke',"black",'stroke-width',1);
			}

			$svg->text('x',$width/2,'y',$height-10,'font-family','TimesNewRoman','font-size',30,'font-color',"black",'text-anchor','middle')->cdata("$query($strand)");
			$svg->rect('x',$margin,'y',$height-$margin,'width',$width-2*$margin,'height',20,'fill',"white",'stroke',"black",'stroke-width',1);
			foreach my $pos (keys %qgene) {
				my ($start,$end)=split /\t/,$pos;
				my $color=$qgene{$pos}[0];
				my $name=$qgene{$pos}[1];
				$svg->text('x',$margin+$start*$xpels2+($end-$start+1)/2*$xpels2,'y',$height-$margin+38,'font-family','TimesNewRoman','font-size',10,'font-color',"$color",'writing-mode','tb','text-anchor','start')->cdata("$name");
				$svg->rect('x',$margin+$start*$xpels2,'y',$height-$margin,'width',($end-$start+1)*$xpels2,'height',20,'fill',"$color",'stroke',"$color",'stroke-width',1);
			}
			$svg->text('x',$margin,'y',$height-$margin+35,'font-family','TimesNewRoman','font-size',10,'font-color',"black",'text-anchor','end')->cdata("0");
			$svg->line('x1',$margin,'y1',$height-$margin+20,'x2',$margin,'y2',$height-$margin+25,'stroke',"black",'stroke-width',1);
			$svg->text('x',$width-$margin,'y',$height-$margin+35,'font-family','TimesNewRoman','font-size',10,'font-color',"black",'text-anchor','start')->cdata( (sprintf "%.2e",$query_len) );
			$svg->line('x1',$width-$margin,'y1',$height-$margin+20,'x2',$width-$margin,'y2',$height-$margin+25,'stroke',"black",'stroke-width',1);
			for (my $i=1;$i<$query_len/20;$i++) {
				$svg->text('x',$margin+$i*($width-2*$margin)/20,'y',$height-$margin+35,'font-family','TimesNewRoman','font-size',10,'font-color',"black",'text-anchor','middle')->cdata( (sprintf "%.2e",$i*($width-2*$margin)/20/$xpels2) );
				$svg->line('x1',$margin+$i*($width-2*$margin)/20,'y1',$height-$margin+20,'x2',$margin+$i*($width-2*$margin)/20,'y2',$height-$margin+25,'stroke',"black",'stroke-width',1);
			}
		}
		my $color;
		if (abs($list[2]-$list[1])==abs($list[4]-$list[3])) {
			$color="blue";
		}else {
			$color="purple";
		}
		my $str="M".($margin+$list[1]*$xpels1)." ".($margin)." L".($margin+$list[2]*$xpels1)." ".($margin)." L".($margin+$list[4]*$xpels2)." ".($height-$margin)." L".($margin+$list[3]*$xpels2)." ".($height-$margin)." Z";
		$svg->path('d',$str,'fill',"$color",'stroke',"$color",'stroke-width',0);
	}
	print O $svg->xmlify();
	close O;
	close IN;
	return;
}

sub tbl2exon {
	my $infile=shift;
	
	my %gene=();
	my $target;
	my $name="";
	my $gene_start="";
	my $gene_end="";
	my $exon_num=0;
	open (IN,"<$infile")||die("fail to open $infile.\n");
	while (<IN>) {
		chomp;
		if (/^>\S+\s+(\S+)$/) {
			$target=$1;
		}elsif (/^(\d+)\t(\d+)\tgene$/) {
			$name="";$exon_num=0;
			if ($1<$2) {
				$gene_start=$1;
				$gene_end=$2;
			}else {
				$gene_start=$2;
				$gene_end=$1;
			}
		}elsif (/^\t\t\tgene\t(.*)$/) {
			$name=$1;
		}elsif (/^(\d+)\t(\d+)\tCDS$/) {
			my $pos1=$1;
			my $pos2=$2;
			$exon_num++;
			if ($pos1<$pos2) {
				$gene{"$pos1\t$pos2"}[0]="red";
				$gene{"$pos1\t$pos2"}[1]=$name;
			}else {
				$gene{"$pos2\t$pos1"}[0]="green";
				$gene{"$pos2\t$pos1"}[1]=$name;
			}
		}elsif (/^(\d+)\t(\d+)\ttRNA$/) {
			$name="tRNA";
			my $pos1=$1;
			my $pos2=$2;
			$exon_num++;
			if ($pos1<$pos2) {
				$gene{"$pos1\t$pos2"}[0]="red";
				$gene{"$pos1\t$pos2"}[1]=$name;
			}else {
				$gene{"$pos2\t$pos1"}[0]="green";
				$gene{"$pos2\t$pos1"}[1]=$name;
			}
		}elsif (/^(\d+)\t(\d+)\trRNA$/) {
			my $pos1=$1;
			my $pos2=$2;
			$name="rRNA";
			$exon_num++;
			if ($pos1<$pos2) {
				$gene{"$pos1\t$pos2"}[0]="red";
				$gene{"$pos1\t$pos2"}[1]=$name;
			}else {
				$gene{"$pos2\t$pos1"}[0]="green";
				$gene{"$pos2\t$pos1"}[1]=$name;
			}
		}elsif (/^(\d+)\t(\d+)$/) {
			my $pos1=$1;
			my $pos2=$2;
			$exon_num++;
			if ($pos1<$pos2) {
				$gene{"$pos1\t$pos2"}[0]="red";
				$gene{"$pos1\t$pos2"}[1]=$name;
			}else {
				$gene{"$pos2\t$pos1"}[0]="green";
				$gene{"$pos2\t$pos1"}[1]=$name;
			}
		}elsif (/^\t\t\tproduct\t(.*)$/) {
			next;
		}
	}
	close IN;
	return \%gene;
}

sub gff2exon {
	my $infile=shift;
	
	my %gene=();
	open (IN,"<$infile")||die("fail to open $infile.\n");
	my $trannew="";
	while (<IN>) {
		chomp;
		#chr1	Cufflinks	exon	4886744	4886831	.	+	.	gene_id "XLOC_000005"; transcript_id "TCONS_00000014"; exon_number "4"; oId "CUFF.13.9"; tss_id "TSS8";
		if (/^\#/) {
			next;
		}
		my @list=split/\t/,$_;
		if ($list[2] ne "exon") {
			next;
		}		
		if ($list[8]=~/transcript_id=([^\;]+)/) {
			$trannew=$1;
		}elsif ($list[8]=~/transcript_id\s+\"([^\"]+)/) {
			$trannew=$1;
		}
		if ($list[6] eq "+") {
			$gene{"$list[3]\t$list[4]"}[0]="red";
		}else {
			$gene{"$list[3]\t$list[4]"}[0]="green";
		}
		$gene{"$list[3]\t$list[4]"}[1]=$trannew;
	}
	close IN;
	return \%gene;
}

sub bed2exon {
	my $infile=shift;

	my %gene=();
	open (IN,"<$infile")||die("fail to open $infile.\n");
	while(<IN>) {
		chomp;
		if (/^track|^\#/) {
			next;
		}
		#chr1	66999065	67210057	ENST00000237247	0	+	67000041	67208778	0	27	25,123,64,25,84,57,55,176,12,12,25,52,86,93,75,501,81,128,127,60,112,156,133,203,65,165,1302,	0,863,92464,99687,100697,106394,109427,110161,127130,134147,137612,138561,139898,143621,146295,148486,150724,155765,156807,162051,185911,195881,200365,205952,207275,207889,209690,
		my @list=split /\t/,$_;
		my @len=split /\,/,$list[10];
		my @start=split /\,/,$list[11];
		for (my $i=0;$i<@start;$i++) {
			if ($list[5] eq "+") {
				$gene{($list[1]+$start[$i]+1)."\t".($list[1]+$start[$i]+$len[$i])}[0]="red";
			}else {
				$gene{($list[1]+$start[$i]+1)."\t".($list[1]+$start[$i]+$len[$i])}[0]="green";
			}
			$gene{($list[1]+$start[$i]+1)."\t".($list[1]+$start[$i]+$len[$i])}[1]=$list[3];
		}
	}
	close IN;
	return \%gene;
}

