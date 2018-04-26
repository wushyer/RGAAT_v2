#!/usr/bin/perl
#Copyright (c) BIG
#Author:Liu Wanfei <liuwf@big.ac.cn>
#Date:2015Äê11ÔÂ02ÈÕ
#Description:This program can filter the variant according to the vcf file.
#use strict;
#use warnings;
my $version="1.0 version";
use Getopt::Long;
#use List::Util qw[min max];
my %opts;
GetOptions(\%opts,"i:s","vc:s","rc:s","o:s");
print "*************\n*$version*\n*************\n";
if (!defined $opts{i} || !defined $opts{o}) {
	die "************************************************************************
	Usage: variant_titv_filter.pl -i vcf_file -o outfile.
	   -i: the absolute path of vcf file.
	  -vc: coverage cutoff of variant.
	  -rc: coverage cutoff of read.
	   -o: the absolute path of outfile.
************************************************************************\n";
}
my $vcffile=$opts{i};
my $outfile=$opts{o};

my (@variant,@ti,@rdp,@dp,@cv,@cr,@bq,@bqb,@bqa,@pc,@mc,@sc,@ic,@dc,@gq,@pl,@mp,);

my $cov_var=1;
if (defined $opts{vc}) {
	$cov_var=$opts{vc};
}
my $cov_read=1;
if (defined $opts{rc}) {
	$cov_read=$opts{rc};
}
my %Ti = (
	"CT" => "0",
	"TC" => "1",
	"AG" => "2",
	"GA" => "3",
);

$time=localtime;
print "Start time: $time.\n";

my $all=0;
my $filtered=0;
open (OUT,">$outfile")||die("fail to open $outfile.\n");
open (IN,"<$vcffile")||die("fail to open $vcffile.\n");
while (<IN>) {
	if (/^\#/) {
		print OUT "$_";
		next;
	}
	$all++;
	my @list=split /\t/,$_;
	my %attr=();
	$attr{"REF"}=$list[3];
	$attr{"ALT"}=$list[4];
	my @ALT=split /\,/,$list[4];
	for (my $i=0;$i<@ALT;$i++) {
		$attr{"ALT$i"}=$ALT[$i];
	}
	$attr{"QUAL"}=$list[5];
	$attr{"FILTER"}=$list[6];
	my @attributes = split /;/, $list[7];
	foreach my $attr ( @attributes) {
		if ($attr =~ /^(\S+)\=(\S+)$/) {
			my $c_type  = $1;
			my $c_value = $2;
			$attr{$c_type} = $c_value;
			my @c_value=split /\,/,$c_value;
			for (my $i=0;$i<@c_value;$i++) {
				$attr{"$c_type$i"}=$c_value[$i];
			}
		}elsif ($attr =~ /^(\S+)$/) {
			my $c_type  = $1;
			$attr{$c_type} = "";
		}
	}
	my @attributes1 = split /:/, $list[8];
	my @attributes2 = split /:/, $list[9];
	for (my $i=0;$i<@attributes1;$i++) {
		my $c_type  = $attributes1[$i];
		my $c_value = $attributes2[$i];
		$attr{$c_type} = $c_value;
		my @c_value=split /\,/,$c_value;
		for (my $j=0;$j<@c_value;$j++) {
			$attr{"$c_type$j"}=$c_value[$j];
		}
	}
	#next if ($attr{"BQ"}<20 and $attr{"IC"}==0 and $attr{"DC"}==0);
	next if ($attr{"BQ"}<20);
	next if ($attr{"AD0"}+$attr{"AD1"}>$attr{"DP"});
	next if ($attr{"AF"}<0.15);
	next if ($attr{"AD0"}/$attr{"DP"}>=0.80);
	if ($attr{"MP"}>0) {
		next if (($attr{"AD1"}-$attr{"MP"})/$attr{"DP"}<0.15);
	}
	if ($attr{"DP"}==1) {
		next;
	}elsif ($attr{"DP"}>1) {
		next if ($attr{"SC"}<2);
	}
	if ($attr{"AD1"}<=$attr{"AD0"}) {
		if ($attr{"AC"}==1) {
			my @min=sort {$a <=> $b} ($attr{"PL0"},$attr{"PL2"});
			next if ($min[0]<=$attr{"PL1"});
		}elsif ($attr{"AC"}==2) {
			my @min=sort {$a <=> $b} ($attr{"PL0"},$attr{"PL1"});
			next if ($min[0]<=$attr{"PL2"});
		}
	}
	$filtered++;
	push (@variant,$_);
	push (@rdp,$attr{"RDP"});
	push (@dp,$attr{"DP"});
	push (@cr,$attr{"AD0"});
	push (@cv,$attr{"AD1"});
	push (@bq,$attr{"BQ"});
	push (@bqb,$attr{"BQB"});
	push (@bqa,$attr{"BQA"});
	push (@pc,$attr{"PC"});
	push (@mc,$attr{"MC"});
	push (@sc,$attr{"SC"});
	push (@ic,$attr{"IC"});
	push (@dc,$attr{"DC"});
	push (@mp,$attr{"MP"});
	if (exists $Ti{"$attr{\"REF\"}$attr{\"ALT0\"}"}) {
		push (@ti,1);
	}else {
		push (@ti,0);
	}
	if ($attr{"AC"}==1) {
		my @min=sort {$a <=> $b} ($attr{"PL0"},$attr{"PL2"});
		push (@gq,$min[0]);
		push (@pl,$attr{"PL1"});
	}elsif ($attr{"AC"}==2) {
		my @min=sort {$a <=> $b} ($attr{"PL0"},$attr{"PL1"});
		push (@gq,$min[0]);
		push (@pl,$attr{"PL2"});
	}
}
close IN;

print "All variant is $all.\n";
print "First filtered variant is $filtered.\n";

$filtered=0;
#filter false positive variant
my $window=1000;
for ($i=0;$i<@ti;$i+=$window) {
	for (my $j=$i;$j<$i+$window;$j++) {
		if ($cv[$j]>=3) {
			if ($ic[$j]>=$dp[$j]*0.50 or $dc[$j]>=$dp[$j]*0.50) {
				$filtered++;
				print OUT "$variant[$j]";
				next;
			}
			if ($cv[$j]<=$cr[$j]) {
				if ($cv[$j]/($cv[$j]+$cr[$j])>0.2) {
					$filtered++;
					print OUT "$variant[$j]";
				}elsif ($bq[$j]>=30 and $bq[$j]>=($bqb[$j]+$bqa[$j])/2) {
					$filtered++;
					print OUT "$variant[$j]";
				}
			}else {
				$filtered++;
				print OUT "$variant[$j]";
			}
		}elsif ($cv[$j]==2) {
			if ($ic[$j]>=$dp[$j]*0.50 or $dc[$j]>=$dp[$j]*0.50) {
				if ($bq[$j]>=30) {
					$filtered++;
					print OUT "$variant[$j]";
					next;
				}
			}
			if ($dp[$j]==2) {
				if ($rdp[$j]==2) {
					#if ($bq[$j]>=35 and $bq[$j]>=($bqb[$j]+$bqa[$j])/2) {
					if ($bq[$j]>=35) {
						$filtered++;
						print OUT "$variant[$j]";
					}
				}elsif ($rdp[$j]>=3) {
					#if ($bq[$j]>=30 and $bq[$j]>=($bqb[$j]+$bqa[$j])/2) {
					if ($bq[$j]>=30) {
						$filtered++;
						print OUT "$variant[$j]";
					}
				}
			}elsif ($dp[$j]>=3) {
				if ($cv[$j]<=$cr[$j]) {
					if ($cv[$j]/($cv[$j]+$cr[$j])>=0.25) {
						#if ($bq[$j]>=30 and $bq[$j]>=($bqb[$j]+$bqa[$j])/2) {
						if ($bq[$j]>=30) {
							$filtered++;
							print OUT "$variant[$j]";
						}
					}elsif ($bq[$j]>=35 and $bq[$j]>=($bqb[$j]+$bqa[$j])/2) {
						$filtered++;
						print OUT "$variant[$j]";
					}
				}else {
					if ($bq[$j]>=30) {
						$filtered++;
						print OUT "$variant[$j]";
					}
				}
			}
		}
	}
}
close OUT;
print "Final filtered variant is $filtered.\n";

$time=localtime;
print "End time: $time.\n";
exit;

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

sub fisher {
	my $n11=shift;
	my $n12=shift;
	my $n21=shift;
	my $n22=shift;

	my $np1=$n11+$n21;
	my $np2=$n12+$n22;
	my $n1p=$n11+$n12;
	my $n2p=$n21+$n22;
	my $npp=$np1+$np2;
	my $fnp1=&factorial($np1);
	my $fnp2=&factorial($np2);
	my $fn1p=&factorial($n1p);
	my $fn2p=&factorial($n2p);
	my $fnpp=&factorial($npp);
	my $fn11=&factorial($n11);
	my $fn12=&factorial($n12);
	my $fn21=&factorial($n21);
	my $fn22=&factorial($n22);
	my $p=($fnp1*$fnp2*$fn1p*$fn2p)/($fnpp*$fn11*$fn12*$fn21*$fn22)*2;
	#print "$p\n";
	$p=log10($p);
	#print "$p\n";
	#$p=($p)*(-10);
	#print "$p\n";
	$p=(sprintf "%.2f",$p);
	#print "$p\n";
	return ($p);
}

sub factorial{
	my $n=shift;
	if ($n == 0) {
		return 1;
	}elsif ($n == 1) {
		return 1;
	}else {
		return $n * factorial($n - 1);
	}
}

sub log10 {
	my $n = shift;
	return log($n)/log(10);
}
