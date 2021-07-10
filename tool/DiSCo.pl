#!/usr/bin/perl

use strict;
use warnings;
use File::Basename;
use Getopt::Std;
use Cwd qw(getcwd);

my $date = localtime();
my $path = getcwd;
my $command = "$0 @ARGV";



my %opts;
getopts('vhi:p:s:n:d:oqm:f:t:', \%opts);
my $program = basename($0);


if($opts{v}){
	print "$program version 1.0.0\n\n";
	exit(0);
}
if($opts{h}){
	usage();
	exit(0);
}
if(!$opts{i}){
	warn "\nERROR:\tSequence input file missing!\n\n\n";
	usage();
	exit(1);
}
print_all("Start\t$date");
print_all($command);

my $seq = $opts{i};
my $out_p;
if($opts{p}){
	$out_p = $opts{p};
}else{
	$out_p = fileparse($seq,'\..*');
}



my $hmm = "$path/DiSCo_HMMlib";
if($opts{m}){
	$hmm=$opts{m};
}
if(!-e $hmm){
	die "Couldn't find $hmm\n";
}


#number of threads to run hmmsearch
my $n = 1;
if($opts{n}){
	$n = $opts{n};
}


#output directory
my $path_o;
if($opts{d}){
	$path_o = $opts{d};
	if($path_o=~m/^(.+)\/$/){
		$path_o=$1;
	}
}else{
	$path_o=$path;	
}


my $hmmout = "$path_o/$out_p.DiSCo.txt";



my $hmmsearch = "hmmsearch -o /dev/null --cpu $n --domtblout $hmmout $hmm $seq ";
print_all($hmmsearch);
my $status = system($hmmsearch);
if($status >0){
	warn "ERROR: $hmmsearch\nDied with exit code:\t$status\n";
	exit($status);
}
my $out = "$path_o/$out_p.DiSCo.filtered.txt";

my $filter = "$path/filter_DiSCo.pl";
if($opts{f}){
	$filter=$opts{f};
}

if(!-e $filter){
	die "Couldn't find $filter\n";
}



my $cmd = "perl $filter -f $hmmout -s $seq -t $out";

#hitlist format
if($opts{s} && $opts{s}<4){
	my $sep =$opts{s};
	$cmd .=" -d $sep";
}elsif($opts{s} && $opts{s}>=4){
	warn "Option -s $opts{s} is not between 0 and 3\nContinue with defaut separator...\n"
}


if($opts{o}){
	my $seq_out = "$path_o/$out_p.DiSCo.faa";
	$cmd .="  -o $seq_out";
}

if($opts{t}){
	my $cut = $opts{t};
	if(!-e $cut){
		warn "$cut does not exist!\nContinue with default thresholds...\n";
	}else{
		$cmd .="  -c $cut";
		warn "$cut: no enzyme-type prediction possible!\n"
	}
}


print_all($cmd);
$status = system($cmd);
if($status >0){
	warn "ERROR: $filter\nDied with exit code:\t$status\n";
	exit($status);
}


#system("mv $out $hmmout");


$date = localtime();
print_all("Done\t$date");



### SUBROUTINES ###
###################


### usage ###
sub usage {
print <<END_OF_USAGE;

DESCRIPTION:
   DiSCo: Dsr-dependent dIssimilatory Sulfur metabolism Classification tOol

USAGE:
  ./$program -i <filename> \[-p <file_prefix>\] \[-o\]
  \[-d </path/to/dir/>\]
  \[-m </path/filename>\] \[-f </path/filename>\]
  \[-t </path/filename>\]
  \[-s <number>\] \[-n <number>\] 
  \[-q\] \[-h\] \[-v\]

 example:
  perl $program -i example.faa  (runs DiSCo with default settings)

  perl $program -i example.faa -p example  (runs DiSCo with given filename prefix)
  perl $program -i example.faa -p example -s 1 -o (runs DiSCo with given filename prefix, specifies outfile separator, and produces sequence outfile)

Standard settings:
  -i <filename>\t\t\tinput protein sequence file in FASTA format


Additional settings:
  -p <prefix>\t\t\tfiltered DiSCo hitlist - filename prefix
  -s <number>\t\t\tseparator for hitlist
   0 = TAB=\t\t\t[default: tab separated]
   1 = ,
   2 = ;
   3 = SPACE=" "
	
  -d </path/to/dir>\t\tabsolute path to outfile directory\t\[optional, default current dir:$path\]
  -o \t\t\t\tDiSCo hits - sequences output in FASTA format\t\[optional\]
  -n <integer>\t\t\tnumber of threads for hmmsearch\t\t\[optional, default = 1\]
  -m \t\t\t\talternative location and filename of DiSCo HMM library\t\[optional, default current dir:$path\]
  -f \t\t\t\talternative location and filename of filter script\t\[optional, default current dir:$path\]
  -t </path/thresholds.txt>\ttab seperated table with user-defined cut-offs: column 1: model, column 2: score, column 3: E-value

General settings:
  -h \t\t\t\tshow help and usage
  -q \t\t\t\tquiet, non verbose
  -v \t\t\t\tversion


END_OF_USAGE
}

sub print_all{
	my $string=shift;
	if(!$opts{q}){
		print "$string\n";
	}
}

