#!/home/wusy/miniconda2/bin/perl
#run blat with Multi CPUs
#author: gaoshh@big.ac.cn wushy@big.ac.cn

use File::Spec::Functions qw/rel2abs/;
use File::Basename;
use Getopt::Long;
use threads;
use Bio::Seq;
use Bio::SeqIO;

#get options
my %opts;
GetOptions(\%opts,
	"h",
	"q:s","t:s","o:s",
	"w:s","p:s"
);

if ((!$opts{q})||(!$opts{t})||(!$opts{o})||($opts{h}))
{
	print "
usage: $0 -h -q <query> -t <target> -o output_psl -p num_cpus -- other BLAT params

Main options:
  -h                     show the help info
  -q   <query>           BLAT query
  -t   <target>          BLAT target
  -o   output_psl        the output file  

  -w   work_path         specify the path where the script works
  -p   cpu               how many parts to split(default 4)

  --                     other BLAT params could be passed in with --, see BLAT manual for more info
";
	exit(1);
}
#settings


#initialize parameters
($QUERY,$TARGET,$OUTPUT,$WORKDIR,$CPUS)=@opts{'q','t','o','w','p'};
#wrap BLAT params
$PARAMS=join " ",@ARGV;

#default
$WORKDIR = ($WORKDIR)?$WORKDIR:rel2abs("./");	##current dir as workdir
$CPUS = ($CPUS > 0)? $CPUS : 4;	#use 4 cpus


#check parameters and environment
$QUERY=rel2abs($QUERY);
$TARGET=rel2abs($TARGET);
unless (-e $QUERY) {die "The QUERY file doesn't exist or unreadable: $QUERY\n";}
unless (-e $TARGET) {die "The TARGET file doesn't exist or unreadable: $TARGET\n";}
if (system("which blat")){
	die "Can't find blat in your system\n";
}
$BLAT = `which blat`;
chomp($BLAT);

unless (-e $WORKDIR) {mkdir $WORKDIR;}
$WORKDIR=rel2abs($WORKDIR);
#map{mkdir "$WORKDIR/$_";}('split','worker');

#link file
system("ln -sf $QUERY $WORKDIR/query.fa");
system("ln -sf $TARGET $WORKDIR/target.fa");


#split input fasta
splitter("$WORKDIR/query.fa",$CPUS) || die "Split input fasta error\n";



#generate shell command
my @cmd;
for $i(1..$CPUS){
	$cmd="$BLAT -noHead $WORKDIR/target.fa $WORKDIR/query.fa.split.$i $WORKDIR/output.split.$i $PARAMS > $WORKDIR/log.$i\n";
	push @cmd,$cmd;
}

#wrap
my @threads;
for $i(1..$CPUS){
	$cmd = shift @cmd;
	$thread = threads->create(\&wrapper,$cmd);
	$threads[$i] = $thread;
}

#recycle 
my $failcount;
my $ret;
for $i(1..$CPUS){
	if ($ret=$threads[$i]->join()){
		$failcount++;
		print "Thread $i exited with non-zero $ret, your result may not be complete, suggest retry...\n";
	}
}

#cat output
die "Some parts of the work end up with error, please check and retry!\n" if ($failcount);
$cmd="cat ";
for $i(1..$CPUS){
	$cmd.="$WORKDIR/output.split.$i ";
}
$cmd.="> $OUTPUT";
system $cmd;
print "Blat finished! Your result could be found: $OUTPUT\n";


#thread
sub wrapper{
	#in: cmd
	#out: retcode
	my $cmd = shift @_;
	print "start CMD: $cmd\n";
	return system($cmd);
}


sub splitter{
	my $fasta=shift @_;
	my $parts=shift @_;

	for $n(1..$parts){
		open $n,">","$fasta.split.$n" || die;
	}
	
	open FASTA,"<",$fasta || die;
	my $l=0;
	while(<FASTA>){
		if (/^>/){
			$l++;
			$n=$l%$parts;$n=$parts if ($n==0);
		}
		print $n $_;
	}

	for $n(1..$parts){
		close OUT;
	}
	close FASTA;
}
