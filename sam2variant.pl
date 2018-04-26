#!/usr/bin/perl
#Copyright (c) BIG
#Author:Liu Wanfei <liuwf@big.ac.cn>
#Date:2015Äê11ÔÂ02ÈÕ
#Description:This program can identify variant according to the sam file.
#use strict;
#use warnings;
my $version="1.0 version";
use Getopt::Long;
#use List::Util qw[min max];
my %opts;
GetOptions(\%opts,"s:s","r:s","cf:s","f:s","cv:s","cr:s","o:s");
print "*************\n*$version*\n*************\n";
if (!defined $opts{r} || !defined $opts{o}) {
	die "************************************************************************
	Usage: sam2variant.pl -s sam_file -r reference -cf cutoff_for_allele_frequency(0~1) -f filter_unpaired_reads -cv coverage_cutoff_of_variant -cr coverage_cutoff_of_read -o outfile.
	   -s: the absolute path of sam file. If there is no -s, then program will read file from standard input.
	  -cf: the cutoff value for allele frequence ranging form 0 to 1.
	   -f: filter unpaired reads.
	  -cv: coverage cutoff of variant.
	  -cr: coverage cutoff of read.
	   -r: the absolute path of reference file.
	   -o: the absolute path of outfile.
************************************************************************\n";
}
my $samfile=$opts{s};
my $reffile=$opts{r};
my $outfile=$opts{o};

#mismatch cutoff
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
my %Ti = (
	"CT" => "0",
	"TC" => "1",
	"AG" => "2",
	"GA" => "3",
);
my $freq_cut=0.15;
if (defined $opts{cf}) {
	$freq_cut=$opts{cf};
}
my $cov_var=1;
if (defined $opts{cv}) {
	$cov_var=$opts{cv};
}
my $cov_read=1;
if (defined $opts{cr}) {
	$cov_read=$opts{cr};
}

my (%ref,%var,);

$time=localtime;
print "Start time: $time.\n";

open (OUT,">$outfile")||die("fail to open $outfile.\n");
print OUT "##fileformat=VCFv4.1\n";
print OUT "##QUAL, it is the odds ratio of being a true variant versus being false according to the reciprocal of Chi-square value for association between two categorical variables\n";
print OUT "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
print OUT "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n";
print OUT "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth (reads with MQ=255 or with bad mates are filtered)\">\n";
print OUT "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n";
print OUT "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, Chi-square value for association between two categorical variables for genotypes as defined in the VCF specification\">\n";

print OUT "##INFO=<ID=RDP,Number=1,Type=Integer,Description=\"All mapped reads depth after reads filter\">\n";
print OUT "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"High quality reads depth; some reads and loci may have been filtered\">\n";
print OUT "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency, for each ALT allele, in the same order as listed\">\n";
print OUT "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes, for each ALT allele, in the same order as listed\">\n";
print OUT "##INFO=<ID=FS,Number=1,Type=Float,Description=\"Phred-scaled p-value using Fisher's exact test to detect strand bias\">\n";
print OUT "##INFO=<ID=QD,Number=1,Type=Float,Description=\"Variant Confidence/Quality by Depth according to the reciprocal of Chi-square value for association between two categorical variables\">\n";
print OUT "##INFO=<ID=BQ,Number=1,Type=Integer,Description=\"Locus base Quality\">\n";
print OUT "##INFO=<ID=MP,Number=1,Type=Integer,Description=\"Multiple mapped reads depth\">\n";
print OUT "##INFO=<ID=IC,Number=1,Type=Integer,Description=\"Insertion reads count\">\n";
print OUT "##INFO=<ID=DC,Number=1,Type=Integer,Description=\"Deletion reads count\">\n";
print OUT "##INFO=<ID=PC,Number=1,Type=Integer,Description=\"Plus strand reads count\">\n";
print OUT "##INFO=<ID=MC,Number=1,Type=Integer,Description=\"Minus strand reads count\">\n";
print OUT "##INFO=<ID=SC,Number=1,Type=Integer,Description=\"Startpoint count\">\n";

open (IN,"<$reffile")||die("fail to open $reffile.\n");
my $id=undef;
my $str=undef;
while (<IN>) {
	chomp;
	if (/^>(\S+)/) {
		if (defined $id) {
			print OUT "##contig=<ID=$id,length=".(length($str)).">\n";
			$str=undef;
		}
		$id=$1;
		next;
	}
	$str.=$_;
}
close IN;
print OUT "##contig=<ID=$id,length=".(length($str)).">\n";
$str=undef;
print OUT "##reference=file://$reffile\n";
print OUT "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	-\n";
my $indelf=0;#the frequency of before indel
my $indelp=0;#the position of before indel
my $indell=0;#the length of before indel

if (defined $samfile) {
	open (STDIN,"<$samfile")||die("fail to open $samfile.\n");
}
my $read_num=0;
my $read_filtered=0;
my $filter_unmap=0;
my $filter_pcr=0;
my $filter_quality=0;
my $filter_chr=0;
my $filter_distance=0;
my $filter_mult=0;
my $filter_repeat=0;
my $filter_short=0;
my $filter_qual=0;
my $read_length=0;
my $startpos=0;
my $chr_skip=0;
$id=undef;
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

	#print variant
	if (defined $id and $id ne $list[2]) {
		($indelf,$indelp,$indell)=&vcf_print($indelf,$indelp,$indell);
		$indelf=0;$indelp=0;$indell=0;%var=();%ref=();
	}
	$id=$list[2];

	#read chromosome
	if (!exists $ref{$list[2]}) {
		open (IN,"<$reffile")||die("fail to open $reffile.\n");
		my $flag=0;
		while (<IN>) {
			chomp;
			if (/^>(\S+)/) {
				if ($flag==0) {
					if ($1 ne $list[2]) {
						next;
					}else {
						$flag=1;
						next;
					}
				}else {
					last;
				}
			}
			$ref{$list[2]}.=$_ if ($flag==1);
		}
		close IN;
		if (!exists $ref{$list[2]}) {
			$ref{$list[2]}="";
			$chr_skip=1;
		}else {
			$chr_skip=0;
		}
	}

	next if ($chr_skip==1);

	#filter reads by average base quality
	my $qual=0;
	foreach my $quality (split //,$list[10]) {
		$qual+=ord($quality)-33;
	}
	$qual=$qual/(length($list[10]));
	if ($qual<20) {
		$filter_qual++;
		next;
	}

	#filter read according to the mapping quality marked in flag
	my $bin=unpack("B32",pack("N",$list[1]));
	$bin=sprintf "%011d",$bin;
	my @bin=split //,$bin;
	if ($bin[8] eq "1") {
		$filter_unmap++;
		next;
	}
	if ($bin[0] eq "1") {
		$filter_pcr++;
		next;
	}
	if ($bin[1] eq "1") {
		$filter_quality++;
		next;
	}
	if (($bin[7] eq "0") and ($bin[8] eq "0")) {
		if ($list[6] ne "=") {
			$filter_chr++;
			next;
		}
		if ($bin[9] eq "0" or abs($list[7]-$list[3])>32000) {
			$filter_distance++;
			next;
		}
	}else {
		if ($bin[2] eq "1") {
			$filter_mult++;
			next;
		}
	}
	
	#read_length
	#$read_length=length($list[9]) if ($read_num<=100 and $read_length<length($list[9]));
	
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
		&cov_read($_);
		$filter_quality++;
		next;
	}elsif ($miscut{$nm}>length($list[9])) {
		&cov_read($_);
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
			$matchlen+=$1 if ($2 eq "M" or $2 eq "D");
			$SM=$1 if (!defined $SM and $2 eq "M");
			$EM=$1 if ($2 eq "M");
		}
	}

	#filter repeat region
	#if ((substr($ref{$list[2]},$list[3]-1,$matchlen))=~/[acgtnN]/) {
	#	$filter_repeat++;
	#	next;
	#}

	while ($SM<8) {
		my ($SMnew,$clipv,$clipr,$label,$num)=&headclip (\@label,\@num);
		&cov_M($list[2],$list[3],$SM);
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
		my $clipr=0;
		for (my $i=0;$i<@label;$i++) {
			if ($label[$i] eq "S") {
			}elsif ($label[$i] eq "M") {
				$clipr+=$num[$i];
			}elsif ($label[$i] eq "N") {
				$clipr+=$num[$i];
			}elsif ($label[$i] eq "I") {
			}elsif ($label[$i] eq "D") {
				$clipr+=$num[$i];
			}
		}
		&cov_M($list[2],$list[3]+$clipr-$EM,$EM);
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
				$var{$list[2]}{$list[3]+$startr}{"all"}++;
				my $qua=substr($list[10],$clip+$startv,1);
				$qua=ord($qua)-33;
				if ($qua>=0.85*$qual) {
					my $var=substr($list[9],$clip+$startv,1);
					my $ref=substr($ref{$list[2]},$list[3]+$startr-1,1);
					if ($ref=~/[acgtnN]+/ or $var=~/[nN]+/) {
					}else {
						$var=uc($var);$ref=uc($ref);
						$var{$list[2]}{$list[3]+$startr}{"total"}++;
						$var{$list[2]}{$list[3]+$startr}{"qual"}+=$qua;
						$var{$list[2]}{$list[3]+$startr}{"start"}++ if ($list[3] ne $startpos);
						if ($bin[2] eq "1") {
							$var{$list[2]}{$list[3]+$startr}{"mult"}++;
						}
						if ($bin[6] eq "1") {
							$var{$list[2]}{$list[3]+$startr}{"minus"}++;
						}else {
							$var{$list[2]}{$list[3]+$startr}{"plus"}++;
						}
						if ($ref eq $var) {
							$var{$list[2]}{$list[3]+$startr}{"ref"}++;
						}else {
							$var{$list[2]}{$list[3]+$startr}{"$ref\t$var"}++;
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
			foreach my $quality (split //,$qua_str) {
				$qua+=ord($quality)-33;
			}
			$qua=$qua/($num[$i]);
			if ($qua>=0.85*$qual) {
				my $var=substr($list[9],$clip+$startv-1,$num[$i]+1);
				my $ref=substr($ref{$list[2]},$list[3]+$startr-2,1);
				if ($ref=~/[acgtnN]+/ or $var=~/[nN]+/) {
				}else {
					$var=uc($var);$ref=uc($ref);
					$var{$list[2]}{$list[3]+$startr-1}{"$ref\t$var"}++;
					$var{$list[2]}{$list[3]+$startr-1}{"insert"}++;
				}
			}
			$startv+=$num[$i];
		}elsif ($label[$i] eq "D") {
			my $qua_str=substr($list[10],$clip+$startv-1,2);
			my $qua=0;
			foreach my $quality (split //,$qua_str) {
				$qua+=ord($quality)-33;
			}
			$qua=$qua/2;
			if ($qua>=0.85*$qual) {
				my $var=substr($list[9],$clip+$startv-1,1);
				my $ref=substr($ref{$list[2]},$list[3]+$startr-2,$num[$i]+1);
				if ($ref=~/[acgtnN]+/ or $var=~/[nN]+/) {
				}else {
					$var=uc($var);$ref=uc($ref);
					$var{$list[2]}{$list[3]+$startr-1}{"$ref\t$var"}++;
					$var{$list[2]}{$list[3]+$startr-1}{"delete"}++;
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
if (defined $samfile) {
	close STDIN;
}

#last chromosome
($indelf,$indelp,$indell)=&vcf_print($indelf,$indelp,$indell);
close OUT;

print "#Total reads is $read_num.\n#Total filtered reads is $read_filtered\n";
print "#The filter reads are as following:\n";
print "Unmap\t$filter_unmap\n";
print "PCR_duplicate\t$filter_pcr\n";
print "Quality_control\t$filter_quality\n";
print "Chromosome_different\t$filter_chr\n";
print "Large_distance\t$filter_distance\n";
print "Multiple\t$filter_mult\n";
print "Repeat_region\t$filter_repeat\n";
print "Possible_error\t$filter_short\n";
print "Low_quality\t$filter_qual\n";

$time=localtime;
print "End time: $time.\n";

sub vcf_print {
	my $indelf=shift;
	my $indelp=shift;
	my $indell=shift;

	foreach my $chr (keys %var) {
		foreach my $pos (sort {$a <=> $b} keys %{$var{$chr}}) {
			my $all=$var{$chr}{$pos}{"all"};
			my $total=$var{$chr}{$pos}{"total"};
			my $qual=0;
			if (defined $var{$chr}{$pos}{"qual"}) {
				$qual=(sprintf "%.0f",$var{$chr}{$pos}{"qual"}/$total) if ($total>0);
			}
			my $refcount=0;
			if (defined $var{$chr}{$pos}{"ref"}) {
				$refcount=$var{$chr}{$pos}{"ref"};
			}
			my $multcount=0;
			if (defined $var{$chr}{$pos}{"mult"}) {
				$multcount=$var{$chr}{$pos}{"mult"};
			}
			my $insertcount=0;
			if (defined $var{$chr}{$pos}{"insert"}) {
				$insertcount=$var{$chr}{$pos}{"insert"};
			}
			my $deletecount=0;
			if (defined $var{$chr}{$pos}{"delete"}) {
				$deletecount=$var{$chr}{$pos}{"delete"};
			}
			my $pluscount=0;
			if (defined $var{$chr}{$pos}{"plus"}) {
				$pluscount=$var{$chr}{$pos}{"plus"};
			}
			my $minuscount=0;
			if (defined $var{$chr}{$pos}{"minus"}) {
				$minuscount=$var{$chr}{$pos}{"minus"};
			}
			my $startcount=0;
			if (exists $var{$chr}{$pos}{"start"}) {
				$startcount=$var{$chr}{$pos}{"start"};
			}
			if ($total<$cov_read) {
				next;
			}
			my $flags=0;
			my $flagi=0;
			foreach my $var (sort { $var{$chr}{$pos}{$b}<=>$var{$chr}{$pos}{$a} } keys %{$var{$chr}{$pos}}) {
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
				next if ($var{$chr}{$pos}{$var}<$cov_var);
				my $varcount=$var{$chr}{$pos}{$var};
				my $AF=(sprintf "%.2f",$varcount/$total);
				next if ($AF<$freq_cut);
				my $AC=2;
				my $GT="1/1";
				my $AD="$refcount,$varcount";

				if (length($ref)==length($read) and $flags==0 and $flagi==0) {
					if ($refcount>0) {
						if ($refcount/$total<0.15) {
							$GT="1/1";
							$AC="2";
						}else {
							$GT="0/1";
							$AC="1";
						}
					}else {
						$GT="1/1";
						$AC="2";
					}
					if ($pos<=$indelp+$indell-1) {
						if ($indelf<$var{$chr}{$pos}{$var}/$total) {
							if ($multcount>0) {
								if (($varcount-$multcount)/$total>=$freq_cut) {
									if (($pluscount>0 and $minuscount>0) or $startcount>=2) {
										my $gt00=& chiX2 (($refcount+$varcount),0,$refcount,$varcount);
										my $gt01=& chiX2 (($refcount+$varcount)/2,($refcount+$varcount)/2,$refcount,$varcount);
										my $gt11=& chiX2 (0,($refcount+$varcount),$refcount,$varcount);
										my $QD=1000;
										my $QUAL=1000;
										my $GQ=0;
										if ($AC==1) {
											my @min=sort {$a <=> $b} ($gt00,$gt11);
											$GQ=$min[0];
											if ($gt01>0) {
												$QD=(sprintf "%.2f",1/$gt01) if (1/$gt01<1000);
												$QUAL=(sprintf "%.2f",$gt00/$gt01) if ($gt00/$gt01<1000);
											}
										}elsif ($AC==2) {
											my @min=sort {$a <=> $b} ($gt00,$gt01);
											$GQ=$min[0];
											if ($gt11>0) {
												$QD=(sprintf "%.2f",1/$gt11) if (1/$gt11<1000);
												$QUAL=(sprintf "%.2f",$gt00/$gt11) if ($gt00/$gt11<1000);
											}
										}
										printf OUT "$chr\t$pos\t\.\t$ref\t$read\t$QUAL\t\.\tRDP=$all;DP=$total;AF=$AF;AC=$AC;QD=$QD;BQ=$qual;MP=$multcount;IC=$insertcount;DC=$deletecount;PC=$pluscount;MC=$minuscount;SC=$startcount\tGT:AD:DP:GQ:PL\t$GT:$AD:$total:$GQ:$gt00,$gt01,$gt11\n";
									}
								}
							}else {
								if (($pluscount>0 and $minuscount>0) or $startcount>=2) {
									my $gt00=& chiX2 (($refcount+$varcount),0,$refcount,$varcount);
									my $gt01=& chiX2 (($refcount+$varcount)/2,($refcount+$varcount)/2,$refcount,$varcount);
									my $gt11=& chiX2 (0,($refcount+$varcount),$refcount,$varcount);
									my $QD=1000;
									my $QUAL=1000;
									my $GQ=0;
									if ($AC==1) {
										my @min=sort {$a <=> $b} ($gt00,$gt11);
										$GQ=$min[0];
										if ($gt01>0) {
											$QD=(sprintf "%.2f",1/$gt01) if (1/$gt01<1000);
											$QUAL=(sprintf "%.2f",$gt00/$gt01) if ($gt00/$gt01<1000);
										}
									}elsif ($AC==2) {
										my @min=sort {$a <=> $b} ($gt00,$gt01);
										$GQ=$min[0];
										if ($gt11>0) {
											$QD=(sprintf "%.2f",1/$gt11) if (1/$gt11<1000);
											$QUAL=(sprintf "%.2f",$gt00/$gt11) if ($gt00/$gt11<1000);
										}
									}
									printf OUT "$chr\t$pos\t\.\t$ref\t$read\t$QUAL\t\.\tRDP=$all;DP=$total;AF=$AF;AC=$AC;QD=$QD;BQ=$qual;MP=$multcount;IC=$insertcount;DC=$deletecount;PC=$pluscount;MC=$minuscount;SC=$startcount\tGT:AD:DP:GQ:PL\t$GT:$AD:$total:$GQ:$gt00,$gt01,$gt11\n";
								}
							}
						}
					}else {
						if ($multcount>0) {
							if (($varcount-$multcount)/$total>=$freq_cut) {
								if (($pluscount>0 and $minuscount>0) or $startcount>=2) {
									my $gt00=& chiX2 (($refcount+$varcount),0,$refcount,$varcount);
									my $gt01=& chiX2 (($refcount+$varcount)/2,($refcount+$varcount)/2,$refcount,$varcount);
									my $gt11=& chiX2 (0,($refcount+$varcount),$refcount,$varcount);
									my $QD=1000;
									my $QUAL=1000;
									my $GQ=0;
									if ($AC==1) {
										my @min=sort {$a <=> $b} ($gt00,$gt11);
										$GQ=$min[0];
										if ($gt01>0) {
											$QD=(sprintf "%.2f",1/$gt01) if (1/$gt01<1000);
											$QUAL=(sprintf "%.2f",$gt00/$gt01) if ($gt00/$gt01<1000);
										}
									}elsif ($AC==2) {
										my @min=sort {$a <=> $b} ($gt00,$gt01);
										$GQ=$min[0];
										if ($gt11>0) {
											$QD=(sprintf "%.2f",1/$gt11) if (1/$gt11<1000);
											$QUAL=(sprintf "%.2f",$gt00/$gt11) if ($gt00/$gt11<1000);
										}
									}
									printf OUT "$chr\t$pos\t\.\t$ref\t$read\t$QUAL\t\.\tRDP=$all;DP=$total;AF=$AF;AC=$AC;QD=$QD;BQ=$qual;MP=$multcount;IC=$insertcount;DC=$deletecount;PC=$pluscount;MC=$minuscount;SC=$startcount\tGT:AD:DP:GQ:PL\t$GT:$AD:$total:$GQ:$gt00,$gt01,$gt11\n";
								}
							}
						}else {
							if (($pluscount>0 and $minuscount>0) or $startcount>=2) {
								my $gt00=& chiX2 (($refcount+$varcount),0,$refcount,$varcount);
								my $gt01=& chiX2 (($refcount+$varcount)/2,($refcount+$varcount)/2,$refcount,$varcount);
								my $gt11=& chiX2 (0,($refcount+$varcount),$refcount,$varcount);
								my $QD=1000;
								my $QUAL=1000;
								my $GQ=0;
								if ($AC==1) {
									my @min=sort {$a <=> $b} ($gt00,$gt11);
									$GQ=$min[0];
									if ($gt01>0) {
										$QD=(sprintf "%.2f",1/$gt01) if (1/$gt01<1000);
										$QUAL=(sprintf "%.2f",$gt00/$gt01) if ($gt00/$gt01<1000);
									}
								}elsif ($AC==2) {
									my @min=sort {$a <=> $b} ($gt00,$gt01);
									$GQ=$min[0];
									if ($gt11>0) {
										$QD=(sprintf "%.2f",1/$gt11) if (1/$gt11<1000);
										$QUAL=(sprintf "%.2f",$gt00/$gt11) if ($gt00/$gt11<1000);
									}
								}
								printf OUT "$chr\t$pos\t\.\t$ref\t$read\t$QUAL\t\.\tRDP=$all;DP=$total;AF=$AF;AC=$AC;QD=$QD;BQ=$qual;MP=$multcount;IC=$insertcount;DC=$deletecount;PC=$pluscount;MC=$minuscount;SC=$startcount\tGT:AD:DP:GQ:PL\t$GT:$AD:$total:$GQ:$gt00,$gt01,$gt11\n";
							}
						}
					}
					$flags=1;
				}

				if (length($ref)!=length($read) and $flagi==0) {
					$indelcount=0;#ref count of indel
					if (length($ref)<length($read)) {
						$indelcount=$total-$insertcount;
					}else {
						$indelcount=$total-$deletecount;
					}
					my $flag=0;
					if ($indelcount>0) {
						if ($indelcount/$total<0.15) {
							$GT="1/1";
							$AC="2";
						}else {
							$GT="0/1";
							$AC="1";
						}
					}else {
						$GT="1/1";
						$AC="2";
					}
					if ($multcount>0) {
						if (($varcount-$multcount)/$total>=$freq_cut) {
							if (($pluscount>0 and $minuscount>0) or $startcount>=2) {
								my $gt00=& chiX2 (($indelcount+$varcount),0,$indelcount,$varcount);
								my $gt01=& chiX2 (($indelcount+$varcount)/2,($indelcount+$varcount)/2,$indelcount,$varcount);
								my $gt11=& chiX2 (0,($indelcount+$varcount),$indelcount,$varcount);
								my $QD=1000;
								my $QUAL=1000;
								my $GQ=0;
								if ($AC==1) {
									my @min=sort {$a <=> $b} ($gt00,$gt11);
									$GQ=$min[0];
									if ($gt01>0) {
										$QD=(sprintf "%.2f",1/$gt01) if (1/$gt01<1000);
										$QUAL=(sprintf "%.2f",$gt00/$gt01) if ($gt00/$gt01<1000);
									}
								}elsif ($AC==2) {
									my @min=sort {$a <=> $b} ($gt00,$gt01);
									$GQ=$min[0];
									if ($gt11>0) {
										$QD=(sprintf "%.2f",1/$gt11) if (1/$gt11<1000);
										$QUAL=(sprintf "%.2f",$gt00/$gt11) if ($gt00/$gt11<1000);
									}
								}
								printf OUT "$chr\t$pos\t\.\t$ref\t$read\t$QUAL\t\.\tRDP=$all;DP=$total;AF=$AF;AC=$AC;QD=$QD;BQ=$qual;MP=$multcount;IC=$insertcount;DC=$deletecount;PC=$pluscount;MC=$minuscount;SC=$startcount\tGT:AD:DP:GQ:PL\t$GT:$indelcount,$varcount:$total:$GQ:$gt00,$gt01,$gt11\n";
								$flag=1;
							}
						}
					}else {
						if (($pluscount>0 and $minuscount>0) or $startcount>=2) {
							my $gt00=& chiX2 (($indelcount+$varcount),0,$indelcount,$varcount);
							my $gt01=& chiX2 (($indelcount+$varcount)/2,($indelcount+$varcount)/2,$indelcount,$varcount);
							my $gt11=& chiX2 (0,($indelcount+$varcount),$indelcount,$varcount);
							my $QD=1000;
							my $QUAL=1000;
							my $GQ=0;
							if ($AC==1) {
								my @min=sort {$a <=> $b} ($gt00,$gt11);
								$GQ=$min[0];
								if ($gt01>0) {
									$QD=(sprintf "%.2f",1/$gt01) if (1/$gt01<1000);
									$QUAL=(sprintf "%.2f",$gt00/$gt01) if ($gt00/$gt01<1000);
								}
							}elsif ($AC==2) {
								my @min=sort {$a <=> $b} ($gt00,$gt01);
								$GQ=$min[0];
								if ($gt11>0) {
									$QD=(sprintf "%.2f",1/$gt11) if (1/$gt11<1000);
									$QUAL=(sprintf "%.2f",$gt00/$gt11) if ($gt00/$gt11<1000);
								}
							}
							printf OUT "$chr\t$pos\t\.\t$ref\t$read\t$QUAL\t\.\tRDP=$all;DP=$total;AF=$AF;AC=$AC;QD=$QD;BQ=$qual;MP=$multcount;IC=$insertcount;DC=$deletecount;PC=$pluscount;MC=$minuscount;SC=$startcount\tGT:AD:DP:GQ:PL\t$GT:$indelcount,$varcount:$total:$GQ:$gt00,$gt01,$gt11\n";
							$flag=1;
						}
					}
					$flagi=1;
					if ($flag=1 and length($read)<length($ref)) {
						$indelf=$var{$chr}{$pos}{$var}/$total;
						$indell=length($ref);
						$indelp=$pos;
					}
				}
				last if ($flags==1 and $flagi==1);
			}
		}
	}
	return ($indelf,$indelp,$indell);
}

sub cov_M {
	my $chr=shift;
	my $start=shift;
	my $length=shift;
	for (my $i=$start;$i<$start+$length;$i++) {
		$var{$chr}{$i}{"all"}++;
	}
}

sub cov_read {
	my $read_str=shift;
	
	my @list=split /\t/,$read_str;
	$list[5]=~s/\d+H+//g;
	my @record=($list[5]=~m/(\d+\D+)/g);
	my @num=();
	my @label=();
	foreach my $record (@record) {
		if ($record=~/(\d+)(\D+)/) {
			push (@num,$1);
			push (@label,$2);
		}
	}
	my $startr=0;
	for (my $i=0;$i<@label;$i++) {
		if ($label[$i] eq "S") {
		}elsif ($label[$i] eq "M") {
			for (my $j=0;$j<$num[$i];$j++) {
				$var{$list[2]}{$list[3]+$startr}{"all"}++;
				$startr++;
			}
		}elsif ($label[$i] eq "N") {
			$startr+=$startr+$num[$i];
		}elsif ($label[$i] eq "I") {
		}elsif ($label[$i] eq "D") {
			$startr+=$num[$i];
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

sub chiX2 {
	my $n11=shift;
	my $n12=shift;
	my $n21=shift;
	my $n22=shift;

	my $np1=$n11+$n21;
	my $np2=$n12+$n22;
	my $n1p=$n11+$n12;
	my $n2p=$n21+$n22;
	my $npp=$np1+$np2;
	#my $m11=$np1*$n1p/$npp;
	#my $m12=$np2*$n1p/$npp;
	#my $m21=$np1*$n2p/$npp;
	#my $m22=$np2*$n2p/$npp;
	my $num1=0;
	my $num2=0;
	my $num3=0;
	my $num4=0;
	$num1=($n11-$n1p/$npp*$np1)**2/($n1p/$npp*$np1) if ($npp*$np1>0);
	$num2=($n12-$n1p/$npp*$np2)**2/($n1p/$npp*$np2) if ($npp*$np2>0);
	$num3=($n21-$n2p/$npp*$np1)**2/($n2p/$npp*$np1) if ($npp*$np1>0);
	$num4=($n22-$n2p/$npp*$np2)**2/($n2p/$npp*$np2) if ($npp*$np2>0);
	my $chisquare=($num1+$num2+$num3+$num4); 
	#print "$chisquare\n";
	#my $phiX2=((($n11*$n22)-($n12*$n21))**2/($n1p*$np1*$np2*$n2p));
	#print "$phiX2\n";
	#my $chiX2=2*((($n11-$m11)/$m11)**2+(($n12-$m12)/$m12)**2+(($n21-$m21)/$m21)**2+(($n22-$m22)/$m22)**2);
	#print "$chiX2\n";
	return (sprintf "%.2f",$chisquare);
}

sub vcf_print2 {
	my $indelf=shift;
	my $indelp=shift;
	my $indell=shift;

	my @variant=();
	my @ti=();
	my @cov=();

	foreach my $chr (keys %var) {
		foreach my $pos (sort {$a <=> $b} keys %{$var{$chr}}) {
			my $total=$var{$chr}{$pos}{"total"};
			my $refcount=0;
			if (defined $var{$chr}{$pos}{"ref"}) {
				$refcount=$var{$chr}{$pos}{"ref"};
			}
			my $multcount=0;
			if (defined $var{$chr}{$pos}{"mult"}) {
				$multcount=$var{$chr}{$pos}{"mult"};
			}
			my $insertcount=0;
			if (defined $var{$chr}{$pos}{"insert"}) {
				$insertcount=$var{$chr}{$pos}{"insert"};
			}
			my $deletecount=0;
			if (defined $var{$chr}{$pos}{"delete"}) {
				$deletecount=$var{$chr}{$pos}{"delete"};
			}
			my $pluscount=0;
			if (defined $var{$chr}{$pos}{"plus"}) {
				$pluscount=$var{$chr}{$pos}{"plus"};
			}
			my $minuscount=0;
			if (defined $var{$chr}{$pos}{"minus"}) {
				$minuscount=$var{$chr}{$pos}{"minus"};
			}
			my $startcount=0;
			if (exists $var{$chr}{$pos}{"start"}) {
				$startcount=$var{$chr}{$pos}{"start"};
			}
			if ($total<$cov_read) {
				next;
			}
			my $flags=0;
			my $flagi=0;
			foreach my $var (sort { $var{$chr}{$pos}{$b}<=>$var{$chr}{$pos}{$a} } keys %{$var{$chr}{$pos}}) {
				next if ($var eq "total");
				next if ($var eq "ref");
				next if ($var eq "mult");
				next if ($var eq "insert");
				next if ($var eq "delete");
				next if ($var eq "plus");
				next if ($var eq "minus");
				next if ($var eq "start");
				my ($ref,$read)=split /\t/,$var;
				next if ($var{$chr}{$pos}{$var}<$cov_var);
				my $varcount=$var{$chr}{$pos}{$var};
				my $AF=(sprintf "%.2f",$varcount/$total);
				next if ($AF<$freq_cut);
				my $AC=2;
				my $GT="1/1";
				my $AD="$refcount,$varcount";

				if (length($ref)==length($read) and $flags==0 and $flagi==0) {
					if ($refcount>0) {
						if ($refcount/$total<0.15) {
							$GT="1/1";
							$AC="2";
						}else {
							$GT="0/1";
							$AC="1";
						}
					}else {
						$GT="1/1";
						$AC="2";
					}
					if ($pos<=$indelp+$indell-1) {
						if ($indelf<$var{$chr}{$pos}{$var}/$total) {
							if ($multcount>0) {
								if (($varcount-$multcount)/$total>=$freq_cut) {
									if (($pluscount>0 and $minuscount>0) or $startcount>=2) {
										my $gt00=& chiX2 (($refcount+$varcount),0,$refcount,$varcount);
										my $gt01=& chiX2 (($refcount+$varcount)/2,($refcount+$varcount)/2,$refcount,$varcount);
										my $gt11=& chiX2 (0,($refcount+$varcount),$refcount,$varcount);
										my $QD=1000;
										my $QUAL=1000;
										if ($AC==1) {
											if ($gt01>0) {
												$QD=(sprintf "%.2f",1/$gt01) if (1/$gt01<1000);
												$QUAL=(sprintf "%.2f",$gt00/$gt01) if ($gt00/$gt01<1000);
											}
										}elsif ($AC==2) {
											if ($gt11>0) {
												$QD=(sprintf "%.2f",1/$gt11) if (1/$gt11<1000);
												$QUAL=(sprintf "%.2f",$gt00/$gt11) if ($gt00/$gt11<1000);
											}
										}
										my @min=sort {$a <=> $b} ($gt00,$gt01,$gt11);
										if (exists $Ti{"$ref$read"}) {
											push (@ti,1);
										}else {
											push (@ti,0);
										}
										push (@cov,$varcount);
										push (@variant,"$chr\t$pos\t\.\t$ref\t$read\t$QUAL\t\.\tDP=$total;AF=$AF;AC=$AC;QD=$QD;BQ=$qual;MP=$multcount;IC=$insertcount;DC=$deletecount;PC=$pluscount;MC=$minuscount;SC=$startcount\tGT:AD:DP:GQ:PL\t$GT:$AD:$total:$min[1]:$gt00,$gt01,$gt11\n");
									}
								}
							}else {
								if (($pluscount>0 and $minuscount>0) or $startcount>=2) {
									my $gt00=& chiX2 (($refcount+$varcount),0,$refcount,$varcount);
									my $gt01=& chiX2 (($refcount+$varcount)/2,($refcount+$varcount)/2,$refcount,$varcount);
									my $gt11=& chiX2 (0,($refcount+$varcount),$refcount,$varcount);
									my $QD=1000;
									my $QUAL=1000;
									if ($AC==1) {
										if ($gt01>0) {
											$QD=(sprintf "%.2f",1/$gt01) if (1/$gt01<1000);
											$QUAL=(sprintf "%.2f",$gt00/$gt01) if ($gt00/$gt01<1000);
										}
									}elsif ($AC==2) {
										if ($gt11>0) {
											$QD=(sprintf "%.2f",1/$gt11) if (1/$gt11<1000);
											$QUAL=(sprintf "%.2f",$gt00/$gt11) if ($gt00/$gt11<1000);
										}
									}
									my @min=sort {$a <=> $b} ($gt00,$gt01,$gt11);
									if (exists $Ti{"$ref$read"}) {
										push (@ti,1);
									}else {
										push (@ti,0);
									}
									push (@cov,$varcount);
									push (@variant,"$chr\t$pos\t\.\t$ref\t$read\t$QUAL\t\.\tDP=$total;AF=$AF;AC=$AC;QD=$QD;BQ=$qual;MP=$multcount;IC=$insertcount;DC=$deletecount;PC=$pluscount;MC=$minuscount;SC=$startcount\tGT:AD:DP:GQ:PL\t$GT:$AD:$total:$min[1]:$gt00,$gt01,$gt11\n");
								}
							}
						}
						$flags=1;
					}else {
						if ($multcount>0) {
							if (($varcount-$multcount)/$total>=$freq_cut) {
								if (($pluscount>0 and $minuscount>0) or $startcount>=2) {
									my $gt00=& chiX2 (($refcount+$varcount),0,$refcount,$varcount);
									my $gt01=& chiX2 (($refcount+$varcount)/2,($refcount+$varcount)/2,$refcount,$varcount);
									my $gt11=& chiX2 (0,($refcount+$varcount),$refcount,$varcount);
									my $QD=1000;
									my $QUAL=1000;
									if ($AC==1) {
										if ($gt01>0) {
											$QD=(sprintf "%.2f",1/$gt01) if (1/$gt01<1000);
											$QUAL=(sprintf "%.2f",$gt00/$gt01) if ($gt00/$gt01<1000);
										}
									}elsif ($AC==2) {
										if ($gt11>0) {
											$QD=(sprintf "%.2f",1/$gt11) if (1/$gt11<1000);
											$QUAL=(sprintf "%.2f",$gt00/$gt11) if ($gt00/$gt11<1000);
										}
									}
									my @min=sort {$a <=> $b} ($gt00,$gt01,$gt11);
									if (exists $Ti{"$ref$read"}) {
										push (@ti,1);
									}else {
										push (@ti,0);
									}
									push (@cov,$varcount);
									push (@variant,"$chr\t$pos\t\.\t$ref\t$read\t$QUAL\t\.\tDP=$total;AF=$AF;AC=$AC;QD=$QD;BQ=$qual;MP=$multcount;IC=$insertcount;DC=$deletecount;PC=$pluscount;MC=$minuscount;SC=$startcount\tGT:AD:DP:GQ:PL\t$GT:$AD:$total:$min[1]:$gt00,$gt01,$gt11\n");
								}
							}
						}else {
							if (($pluscount>0 and $minuscount>0) or $startcount>=2) {
								my $gt00=& chiX2 (($refcount+$varcount),0,$refcount,$varcount);
								my $gt01=& chiX2 (($refcount+$varcount)/2,($refcount+$varcount)/2,$refcount,$varcount);
								my $gt11=& chiX2 (0,($refcount+$varcount),$refcount,$varcount);
								my $QD=1000;
								my $QUAL=1000;
								if ($AC==1) {
									if ($gt01>0) {
										$QD=(sprintf "%.2f",1/$gt01) if (1/$gt01<1000);
										$QUAL=(sprintf "%.2f",$gt00/$gt01) if ($gt00/$gt01<1000);
									}
								}elsif ($AC==2) {
									if ($gt11>0) {
										$QD=(sprintf "%.2f",1/$gt11) if (1/$gt11<1000);
										$QUAL=(sprintf "%.2f",$gt00/$gt11) if ($gt00/$gt11<1000);
									}
								}
								my @min=sort {$a <=> $b} ($gt00,$gt01,$gt11);
								if (exists $Ti{"$ref$read"}) {
									push (@ti,1);
								}else {
									push (@ti,0);
								}
								push (@cov,$varcount);
								push (@variant,"$chr\t$pos\t\.\t$ref\t$read\t$QUAL\t\.\tDP=$total;AF=$AF;AC=$AC;QD=$QD;BQ=$qual;MP=$multcount;IC=$insertcount;DC=$deletecount;PC=$pluscount;MC=$minuscount;SC=$startcount\tGT:AD:DP:GQ:PL\t$GT:$AD:$total:$min[1]:$gt00,$gt01,$gt11\n");
							}
						}
					}
					$flags=1;
				}

				if (length($ref)!=length($read) and $flagi==0) {
					$indelcount=0;#ref count of indel
					if (length($ref)<length($read)) {
						$indelcount=$total-$insertcount;
					}else {
						$indelcount=$total-$deletecount;
					}
					if ($indelcount>0) {
						if ($indelcount/$total<0.15) {
							$AC=2;
							$GT="1/1";
						}else {
							$GT="0/1";
							$AC="1";
						}
					}else {
						$AC=2;
						$GT="1/1";
					}
					if ($multcount>0) {
						if (($varcount-$multcount)/$total>=$freq_cut) {
							if (($pluscount>0 and $minuscount>0) or $startcount>=2) {
								my $gt00=& chiX2 (($indelcount+$varcount),0,$indelcount,$varcount);
								my $gt01=& chiX2 (($indelcount+$varcount)/2,($indelcount+$varcount)/2,$indelcount,$varcount);
								my $gt11=& chiX2 (0,($indelcount+$varcount),$indelcount,$varcount);
								my $QD=1000;
								my $QUAL=1000;
								if ($AC==1) {
									if ($gt01>0) {
										$QD=(sprintf "%.2f",1/$gt01) if (1/$gt01<1000);
										$QUAL=(sprintf "%.2f",$gt00/$gt01) if ($gt00/$gt01<1000);
									}
								}elsif ($AC==2) {
									if ($gt11>0) {
										$QD=(sprintf "%.2f",1/$gt11) if (1/$gt11<1000);
										$QUAL=(sprintf "%.2f",$gt00/$gt11) if ($gt00/$gt11<1000);
									}
								}
								my @min=sort {$a <=> $b} ($gt00,$gt01,$gt11);
								if (exists $Ti{"$ref$read"}) {
									push (@ti,1);
								}else {
									push (@ti,0);
								}
								push (@cov,$varcount);
								push (@variant,"$chr\t$pos\t\.\t$ref\t$read\t$QUAL\t\.\tDP=$total;AF=$AF;AC=$AC;QD=$QD;BQ=$qual;MP=$multcount;IC=$insertcount;DC=$deletecount;PC=$pluscount;MC=$minuscount;SC=$startcount\tGT:AD:DP:GQ:PL\t$GT:$indelcount,$varcount:$total:$min[1]:$gt00,$gt01,$gt11\n");
							}
						}
					}else {
						if (($pluscount>0 and $minuscount>0) or $startcount>=2) {
							my $gt00=& chiX2 (($indelcount+$varcount),0,$indelcount,$varcount);
							my $gt01=& chiX2 (($indelcount+$varcount)/2,($indelcount+$varcount)/2,$indelcount,$varcount);
							my $gt11=& chiX2 (0,($indelcount+$varcount),$indelcount,$varcount);
							my $QD=1000;
							my $QUAL=1000;
							if ($AC==1) {
								if ($gt01>0) {
									$QD=(sprintf "%.2f",1/$gt01) if (1/$gt01<1000);
									$QUAL=(sprintf "%.2f",$gt00/$gt01) if ($gt00/$gt01<1000);
								}
							}elsif ($AC==2) {
								if ($gt11>0) {
									$QD=(sprintf "%.2f",1/$gt11) if (1/$gt11<1000);
									$QUAL=(sprintf "%.2f",$gt00/$gt11) if ($gt00/$gt11<1000);
								}
							}
							my @min=sort {$a <=> $b} ($gt00,$gt01,$gt11);
							if (exists $Ti{"$ref$read"}) {
								push (@ti,1);
							}else {
								push (@ti,0);
							}
							push (@cov,$varcount);
							push (@variant,"$chr\t$pos\t\.\t$ref\t$read\t$QUAL\t\.\tDP=$total;AF=$AF;AC=$AC;QD=$QD;BQ=$qual;MP=$multcount;IC=$insertcount;DC=$deletecount;PC=$pluscount;MC=$minuscount;SC=$startcount\tGT:AD:DP:GQ:PL\t$GT:$indelcount,$varcount:$total:$min[1]:$gt00,$gt01,$gt11\n");
						}
					}
					$flagi=1;
					if (length($read)<length($ref)) {
						$indelf=$var{$chr}{$pos}{$var}/$total;
						$indell=length($ref);
						$indelp=$pos;
					}
				}
				last if ($flags==1 and $flagi==1);
			}
		}
	}

	#filter false positive variant
	my $TiTv=2;
	my $ti=0;
	my $tv=0;
	my $titv=0;
	my $fdr=0;
	my $window=1000;
	for ($i=0;$i<@ti;$i+=$window) {
		$ti=0;
		$vi=0;
		for (my $j=$i;$j<$i+$window;$j++) {
			if ($ti[$j]==1) {
				$ti++;
			}else {
				$tv++;
			}
		}
		if ($tv>0) {
			$titv=$ti/$tv;
		}else {
			$titv=$TiTv;
		}
		$fdr=($titv-0.5)/($TiTv-0.5);
		if ($fdr<1) {
			for (my $j=$i;$j<$i+$window;$j++) {
				if ($cov[$j]>=3) {
					print OUT "$variant[$j]";
				}
			}
		}else {
			for (my $j=$i;$j<$i+$window;$j++) {
				print OUT "$variant[$j]";
			}
		}
	}
	return ($indelf,$indelp,$indell);
}
