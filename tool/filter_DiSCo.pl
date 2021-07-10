#!/usr/bin/perl

use strict;
use warnings;
use FileHandle;
use Getopt::Std;

use File::Basename;

my %opts;
getopts('s:d:hf:klnrf:o:t:c:', \%opts);

if($opts{h}){
	usage ();
	exit(0);
}elsif(!$opts{s} || !$opts{t}  || !$opts{f} ){
	usage ();
	exit(1);
}
my $sep = "\t";
if($opts{d}){
	my @s=("\t", ",", ";", " ");
	$sep =$s[$opts{d}];
}

my ($in, $faa, $db, $fh, $fho, $head, $Cut, $DsrC, %ids, %head2seq, $king, $m, $d, $switch, $hit_score, $hit_evalue, $hit_accuracy, $o, $c);
my ($line, $id, $tmp, $prot, $i, $j, $h_s, $h_e, $d_e, $d_s, $e_e, $e_s, $genome, $name, $out, $file, $m_score, $m_bias);
my ($hmm_len, $dom_len, $evalue, $model, $score, $anno, $slen, $qlen, $cov, $bias, $accuracy, $current_domain, $total_domnumber);
my ($dom_score, $dom_bias, $c_evalue, $i_evalue);
my (@list, @hit, @todo, @s, @b, %best, %hit, %type, %pheno, %meta, %tax, %do, %file2genome, %cys, %evalue, %score, %e_part, %s_part, %anno, %prot, %pred, %id2o, %id2c, %dir);

if($opts{c}){
	$Cut=$opts{c};
	Userinput($Cut);
}else{
	Thresholds();
}



$DsrC = 0;

$in = $opts{f};

$fh = FileHandle->new;
open($fh, "< $in") or die $!;

while($line=<$fh>){
	chomp($line);
	if($line =~ m/^#/){ 
		$tmp=$line;
		next;
	}
	@list = split(/\s+/, $line);
	$switch = 0;
	$id = $list[0];
	$model = $list[3];
	$qlen = $list[2];
	$slen = $list[5];
	$evalue = $list[6];	
	$score = $list[7];
	$bias = $list[8];
	$current_domain = $list[9];
	$total_domnumber = $list[10];
	$c_evalue = $list[11];	
	$i_evalue = $list[12];
	$dom_score = $list[13];
	$dom_bias = $list[14];
	$h_s = $list[15];
	$h_e = $list[16];
	$hmm_len=$h_e-$h_s+1;
	$d_s = $list[17];
	$d_e = $list[18];
	$dom_len = $d_e-$d_s+1;
	$e_s = $list[19];	
	$e_e = $list[20];
	$accuracy = $list[21];

	$j = scalar(@list)-1;
	$anno = join("_", @list[22..$j]);
	$anno{$id}=$anno;
			
	@s=split('\.', $score);
	$m_score=length($s[0]);	
	@b=split('\.', $bias);
	$m_bias=length($b[0]);	




	if($score >= $score{$model}  && $evalue <= $evalue{$model} && $c_evalue <=1e-10 && $i_evalue <=1e-10 && $accuracy >=0.5 && $m_score > $m_bias){
		if(!exists $hit{$id}){ 
			$switch = 1;
		}elsif($score>$best{$id}){
			$switch = 1;
		}else{
			#hit lower
		}
	}

	if($switch > 0){
		$best{$id}=$score;
		$hit{$id}="$qlen$sep$slen$sep$score$sep$bias$sep$evalue$sep$accuracy$sep$hmm_len$sep$c_evalue$sep$i_evalue$sep$dom_score$sep$dom_bias$sep$current_domain$sep$total_domnumber$sep$d_s$sep$d_e";
		$type{$id}=$model;
		if($model =~m/DsrC/){
			$DsrC++;
		}
	}
}
if($tmp !~m/^#\s*\[ok\]/){
	die "$in incomplete!\n$tmp\n";		
}
close($fh);
######################################################################################################
if($DsrC>0){
	%head2seq = extract_seq(\%hit);
	%cys = cysteines(\%head2seq, \%type);
	foreach $id (keys %type){
		if($type{$id} =~m/DsrC/ && !exists $cys{$id}){
			delete($hit{$id});
		}
	}
}

######################################################################################################
if(scalar(keys %hit) == 0){
	warn "$opts{f}\tNo significant hits\n";
	undef(%hit);
	undef(%best);
	undef(%type);
	undef(%anno);
	undef(%prot);
	undef(%cys);
	undef(%head2seq);
	$DsrC=0;
	exit;
}
######################################################################################################
if($opts{o}){
	if(scalar(keys %head2seq) == 0){
		%head2seq = extract_seq(\%hit);
	}
	print_seq(\%head2seq, \%hit, \%cys);
}

######################################################################################################
$out = $opts{t};
open($fho, "> $out") || die "No such file $out\n";

print $fho "#ID$sep"."protein$sep"."enzyme_type$sep"."qlen$sep"."slen$sep"."score$sep"."bias$sep"."evalue$sep"."accuracy$sep"."hmm_len$sep"."c_evalue$sep"."i_evalue$sep"."dom_score$sep"."dom_bias$sep"."current_domain$sep"."total_domnumber$sep"."domain_start$sep"."domain_end$sep"."annotation$sep"."notes\n";

foreach $id (sort {substr($type{$a},0,4) cmp substr($type{$b},0,4) } keys %hit){
	$tmp = "$id$sep$type{$id}$sep$dir{$type{$id}}$sep$hit{$id}$sep$anno{$id}";

	print $fho "$tmp";
	if($type{$id} =~m/DsrC/){
		print $fho "$sep$cys{$id}";
	}
	print $fho "\n";

}
close($fho);
undef(%hit);
undef(%best);
undef(%type);
undef(%anno);
undef(%prot);
undef(%cys);
undef(%head2seq);
undef(%ids);




######################################################################################################
#
#		Subroutines
#
######################################################################################################
sub extract_seq{
	$in =$opts{s};
	%ids = %{$_[0]};


	$fh = FileHandle->new;
	open($fh, "< $in") or die $!;
	
	while($line=<$fh>){
		chomp($line);
		if($line =~ m/>(\S+)/){
			$head = $1;
			if(exists $head2seq{$head}){
				die "duplicated sequence name\t$in\t$line\n";
				delete($head2seq{$head});
			}
		}else{
			if(exists $ids{$head}){
				$head2seq{$head}.=$line;
			}
		}
	}
	close($fh);
	return(%head2seq);
}
######################################################################################################
sub cysteines{
	my %seq = %{$_[0]};
	my %prot = %{$_[1]};
	my $motif="AGLPKPTG";
	my %true;
	foreach $id (keys %seq){
		if($prot{$id} =~m/DsrC/ && $seq{$id} =~m/C.{10,10}C.{0,4}$/){
			$true{$id}="DsrC-2Cys";
		}elsif($prot{$id} =~m/DsrC/ && $seq{$id} =~m/(.{8,8})C.{0,4}$/){
			my $match=$&;
			my $c=0;
			for(my $p=0; $p<length($motif); $p++){
				if(substr($motif,$p,1) eq substr($match,$p,1)){
					$c++;
				}
			}
			if($c>=7 && $match !~m/NC.{0,4}/){
				$true{$id}="DsrC-1Cys";
			}
		}
	}
	return(%true);
}
######################################################################################################
sub print_seq{
	my %seq = %{$_[0]};
	my %ids = %{$_[1]};
	my %true =%{$_[2]};
	$out = $opts{o};
	open(OUT, "> $out") || die "No such file $out\n";
	foreach $id (sort keys %ids){
		my $tmp="";
		if(exists $true{$id}){
			$tmp="$sep".$true{$id};
		}
		print OUT ">$id$sep$type{$id}$tmp\n$seq{$id}\n";
	}
	close(OUT);
}
######################################################################################################
sub Thresholds{
	%evalue=(			
		"AprA_arch"	 =>	"1e-10",
		"AprA_ass"	 =>	"1e-10",
		"AprA_Chlorobi"	 =>	"1e-10",
		"AprA_SOB_M"	 =>	"1e-10",
		"AprA_SOB_Q"	 =>	"1e-10",
		"AprA_SRO"	 =>	"1e-10",
		"AprB_arch"	 =>	"1e-10",
		"AprB_ass"	 =>	"1e-10",
		"AprB_Chlorobi"	 =>	"1e-10",
		"AprB_SOB_M"	 =>	"1e-10",
		"AprB_SOB_Q"	 =>	"1e-10",
		"AprB_SRO"	 =>	"1e-10",
		"AprM"	 =>	"1e-10",
		"DsrA_arch"	 =>	"1e-10",
		"DsrA_SDM"	 =>	"1e-10",
		"DsrA_SOB"	 =>	"1e-10",
		"DsrA_SRO"	 =>	"1e-10",
		"DsrB_arch"	 =>	"1e-10",
		"DsrB_SDM"	 =>	"1e-10",
		"DsrB_SOB"	 =>	"1e-10",
		"DsrB_SRO"	 =>	"1e-10",
		"DsrC_arch"	 =>	"1e-10",
		"DsrC_SDM"	 =>	"1e-10",
		"DsrC_SOB"	 =>	"1e-10",
		"DsrC_SRO"	 =>	"1e-10",
		"DsrD1_SDM"	 =>	"1e-10",
		"DsrD1"	 =>	"1e-10",
		"DsrD2"	 =>	"1e-10",
		"DsrD3"	 =>	"1e-10",
		"DsrD4"	 =>	"1e-10",
		"DsrE"	 =>	"1e-10",
		"DsrF"	 =>	"1e-10",
		"DsrH"	 =>	"1e-10",
		"DsrJ_Chlorobi"	 =>	"1e-10",
		"DsrJ_SDM"	 =>	"1e-10",
		"DsrJ_SOB"	 =>	"1e-10",
		"DsrJ_SRO"	 =>	"1e-10",
		"DsrK_Chlorobi"	 =>	"1e-10",
		"DsrK_MK_SRO"	 =>	"1e-10",
		"DsrK_SDM"	 =>	"1e-10",
		"DsrK_SOB"	 =>	"1e-10",
		"DsrK_SRO"	 =>	"1e-10",
		"DsrL_SDM"	 =>	"1e-10",
		"DsrL_SOB"	 =>	"1e-10",
		"DsrM_Chlorobi"	 =>	"1e-10",
		"DsrM_MK_SRO"	 =>	"1e-10",
		"DsrM_SDM"	 =>	"1e-10",
		"DsrM_SOB"	 =>	"1e-10",
		"DsrM_SRO"	 =>	"1e-10",
		"DsrO_Chlorobi"	 =>	"1e-10",
		"DsrO_SDM"	 =>	"1e-10",
		"DsrO_SOB"	 =>	"1e-10",
		"DsrO_SRO"	 =>	"1e-10",
		"DsrP_Chlorobi"	 =>	"1e-10",
		"DsrP_SDM"	 =>	"1e-10",
		"DsrP_SOB"	 =>	"1e-10",
		"DsrP_SRO"	 =>	"1e-10",
		"DsrT_Chlorobi"	 =>	"1e-10",
		"DsrT_SRO"	 =>	"1e-10",
		"HdrA_SOB"	 =>	"1e-10",
		"HdrA_SRO"	 =>	"1e-10",
		"HdrA1_methano"	 =>	"1e-10",
		"HdrA2_methano"	 =>	"1e-10",
		"HdrB_SOB"	 =>	"1e-10",
		"HdrB_SRO"	 =>	"1e-10",
		"HdrB1_methano"	 =>	"1e-10",
		"HdrB2_methano"	 =>	"1e-10",
		"HdrB2_SOB"	 =>	"1e-10",
		"HdrC_SOB"	 =>	"1e-10",
		"HdrC_SRO"	 =>	"1e-10",
		"HdrC1_methano"	 =>	"1e-10",
		"HdrC2_methano"	 =>	"1e-10",
		"HdrC2_SOB"	 =>	"1e-10",
		"QmoA_SOB"	 =>	"1e-10",
		"QmoA_SRO"	 =>	"1e-10",
		"QmoA2_SRO"	 =>	"1e-10",
		"QmoB_SOB"	 =>	"1e-10",
		"QmoB_SRO"	 =>	"1e-10",
		"QmoB2_SRO"	 =>	"1e-10",
		"QmoC_SOB"	 =>	"1e-10",
		"QmoC_SRO"	 =>	"1e-10",
		"Sat_arch"	 =>	"1e-10",
		"Sat_ass"	 =>	"1e-10",
		"Sat_CF"	 =>	"1e-10",
		"Sat_Chlorobi"	 =>	"1e-10",
		"Sat_Delta"	 =>	"1e-10",
		"Sat_NF"	 =>	"1e-10",
		"Sat_small"	 =>	"1e-10",
		"Sat_SOB"	 =>	"1e-10",
		"Sat_SOB2"	 =>	"1e-10",
		"Sat_SRO"	 =>	"1e-10",
	); 			
				
	%score =(			
		"AprA_arch"	 =>	"399.1",
		"AprA_ass"	 =>	"399.1",
		"AprA_Chlorobi"	 =>	"399.1",
		"AprA_SOB_M"	 =>	"399.1",
		"AprA_SOB_Q"	 =>	"399.1",
		"AprA_SRO"	 =>	"753.8",
		"AprB_arch"	 =>	"99.7",
		"AprB_ass"	 =>	"63.0",
		"AprB_Chlorobi"	 =>	"128.4",
		"AprB_SOB_M"	 =>	"160.4",
		"AprB_SOB_Q"	 =>	"128.4",
		"AprB_SRO"	 =>	"128.4",
		"AprM"	 =>	"123.4",
		"DsrA_arch"	 =>	"155.5",
		"DsrA_SDM"	 =>	"155.5",
		"DsrA_SOB"	 =>	"155.5",
		"DsrA_SRO"	 =>	"155.5",
		"DsrB_arch"	 =>	"155.5",
		"DsrB_SDM"	 =>	"155.5",
		"DsrB_SOB"	 =>	"155.5",
		"DsrB_SRO"	 =>	"155.5",
		"DsrC_arch"	 =>	"139.8",
		"DsrC_SDM"	 =>	"139.8",
		"DsrC_SOB"	 =>	"139.8",
		"DsrC_SRO"	 =>	"139.8",
		"DsrD1_SDM"	 =>	"45.6",
		"DsrD1"	 =>	"45.6",
		"DsrD2"	 =>	"45.6",
		"DsrD3"	 =>	"45.6",
		"DsrD4"	 =>	"45.6",
		"DsrE"	 =>	"134.7",
		"DsrF"	 =>	"110.5",
		"DsrH"	 =>	"101.5",
		"DsrJ_Chlorobi"	 =>	"102.3",
		"DsrJ_SDM"	 =>	"102.3",
		"DsrJ_SOB"	 =>	"85.7",
		"DsrJ_SRO"	 =>	"102.3",
		"DsrK_Chlorobi"	 =>	"388.1",
		"DsrK_MK_SRO"	 =>	"388.1",
		"DsrK_SDM"	 =>	"388.1",
		"DsrK_SOB"	 =>	"388.1",
		"DsrK_SRO"	 =>	"388.1",
		"DsrL_SDM"	 =>	"348.7",
		"DsrL_SOB"	 =>	"452.4",
		"DsrM_Chlorobi"	 =>	"213.3",
		"DsrM_MK_SRO"	 =>	"117.5",
		"DsrM_SDM"	 =>	"213.3",
		"DsrM_SOB"	 =>	"143.5",
		"DsrM_SRO"	 =>	"213.3",
		"DsrO_Chlorobi"	 =>	"223.5",
		"DsrO_SDM"	 =>	"223.5",
		"DsrO_SOB"	 =>	"223.5",
		"DsrO_SRO"	 =>	"223.5",
		"DsrP_Chlorobi"	 =>	"375.8",
		"DsrP_SDM"	 =>	"375.8",
		"DsrP_SOB"	 =>	"225.8",
		"DsrP_SRO"	 =>	"375.8",
		"DsrT_Chlorobi"	 =>	"63.5",
		"DsrT_SRO"	 =>	"58.7",
		"HdrA_SOB"	 =>	"183.4",
		"HdrA_SRO"	 =>	"183.4",
		"HdrA1_methano"	 =>	"183.4",
		"HdrA2_methano"	 =>	"183.4",
		"HdrB_SOB"	 =>	"123.4",
		"HdrB_SRO"	 =>	"153.4",
		"HdrB1_methano"	 =>	"123.4",
		"HdrB2_methano"	 =>	"123.4",
		"HdrB2_SOB"	 =>	"123.4",
		"HdrC_SOB"	 =>	"79.1",
		"HdrC_SRO"	 =>	"79.1",
		"HdrC1_methano"	 =>	"79.1",
		"HdrC2_methano"	 =>	"79.1",
		"HdrC2_SOB"	 =>	"79.1",
		"QmoA_SOB"	 =>	"408.3",
		"QmoA_SRO"	 =>	"408.3",
		"QmoA2_SRO"	 =>	"613.8",
		"QmoB_SOB"	 =>	"449.9",
		"QmoB_SRO"	 =>	"613.8",
		"QmoB2_SRO"	 =>	"449.9",
		"QmoC_SOB"	 =>	"142.8",
		"QmoC_SRO"	 =>	"142.8",
		"Sat_arch"	 =>	"282.7",
		"Sat_ass"	 =>	"282.7",
		"Sat_CF"	 =>	"349.5",
		"Sat_Chlorobi"	 =>	"282.7",
		"Sat_Delta"	 =>	"282.7",
		"Sat_NF"	 =>	"282.7",
		"Sat_small"	 =>	"455.4",
		"Sat_SOB"	 =>	"282.7",
		"Sat_SOB2"	 =>	"282.7",
		"Sat_SRO"	 =>	"282.7",
	); 			
				
	%dir =(			
		"AprA_arch"	 =>	"reductive",
		"AprA_ass"	 =>	"unspecific",
		"AprA_Chlorobi"	 =>	"oxidative",
		"AprA_SOB_M"	 =>	"oxidative",
		"AprA_SOB_Q"	 =>	"oxidative",
		"AprA_SRO"	 =>	"reductive",
		"AprB_arch"	 =>	"reductive",
		"AprB_ass"	 =>	"unspecific",
		"AprB_Chlorobi"	 =>	"oxidative",
		"AprB_SOB_M"	 =>	"oxidative",
		"AprB_SOB_Q"	 =>	"oxidative",
		"AprB_SRO"	 =>	"reductive",
		"AprM"	 =>	"oxidative",
		"DsrA_arch"	 =>	"reductive",
		"DsrA_SDM"	 =>	"reductive_dismutation",
		"DsrA_SOB"	 =>	"oxidative",
		"DsrA_SRO"	 =>	"reductive",
		"DsrB_arch"	 =>	"reductive",
		"DsrB_SDM"	 =>	"reductive_dismutation",
		"DsrB_SOB"	 =>	"oxidative",
		"DsrB_SRO"	 =>	"reductive",
		"DsrC_arch"	 =>	"reductive",
		"DsrC_SDM"	 =>	"reductive_dismutation",
		"DsrC_SOB"	 =>	"oxidative",
		"DsrC_SRO"	 =>	"reductive",
		"DsrD1_SDM"	 =>	"reductive_dismutation",
		"DsrD1"	 =>	"reductive",
		"DsrD2"	 =>	"reductive",
		"DsrD3"	 =>	"reductive",
		"DsrD4"	 =>	"reductive",
		"DsrE"	 =>	"oxidative",
		"DsrF"	 =>	"oxidative",
		"DsrH"	 =>	"oxidative",
		"DsrJ_Chlorobi"	 =>	"oxidative",
		"DsrJ_SDM"	 =>	"reductive_dismutation",
		"DsrJ_SOB"	 =>	"oxidative",
		"DsrJ_SRO"	 =>	"reductive",
		"DsrK_Chlorobi"	 =>	"oxidative",
		"DsrK_MK_SRO"	 =>	"reductive",
		"DsrK_SDM"	 =>	"reductive_dismutation",
		"DsrK_SOB"	 =>	"oxidative",
		"DsrK_SRO"	 =>	"reductive",
		"DsrL_SDM"	 =>	"reductive_dismutation",
		"DsrL_SOB"	 =>	"oxidative",
		"DsrM_Chlorobi"	 =>	"oxidative",
		"DsrM_MK_SRO"	 =>	"reductive",
		"DsrM_SDM"	 =>	"reductive_dismutation",
		"DsrM_SOB"	 =>	"oxidative",
		"DsrM_SRO"	 =>	"reductive",
		"DsrO_Chlorobi"	 =>	"oxidative",
		"DsrO_SDM"	 =>	"reductive_dismutation",
		"DsrO_SOB"	 =>	"oxidative",
		"DsrO_SRO"	 =>	"reductive",
		"DsrP_Chlorobi"	 =>	"oxidative",
		"DsrP_SDM"	 =>	"reductive_dismutation",
		"DsrP_SOB"	 =>	"oxidative",
		"DsrP_SRO"	 =>	"reductive",
		"DsrT_Chlorobi"	 =>	"oxidative",
		"DsrT_SRO"	 =>	"reductive",
		"HdrA_SOB"	 =>	"unspecific",
		"HdrA_SRO"	 =>	"unspecific",
		"HdrA1_methano"	 =>	"unspecific",
		"HdrA2_methano"	 =>	"unspecific",
		"HdrB_SOB"	 =>	"unspecific",
		"HdrB_SRO"	 =>	"unspecific",
		"HdrB1_methano"	 =>	"unspecific",
		"HdrB2_methano"	 =>	"unspecific",
		"HdrB2_SOB"	 =>	"unspecific",
		"HdrC_SOB"	 =>	"unspecific",
		"HdrC_SRO"	 =>	"unspecific",
		"HdrC1_methano"	 =>	"unspecific",
		"HdrC2_methano"	 =>	"unspecific",
		"HdrC2_SOB"	 =>	"unspecific",
		"QmoA_SOB"	 =>	"oxidative",
		"QmoA_SRO"	 =>	"reductive",
		"QmoA2_SRO"	 =>	"reductive",
		"QmoB_SOB"	 =>	"oxidative",
		"QmoB_SRO"	 =>	"reductive",
		"QmoB2_SRO"	 =>	"reductive",
		"QmoC_SOB"	 =>	"oxidative",
		"QmoC_SRO"	 =>	"reductive",
		"Sat_arch"	 =>	"reductive",
		"Sat_ass"	 =>	"unspecific",
		"Sat_CF"	 =>	"unspecific",
		"Sat_Chlorobi"	 =>	"oxidative",
		"Sat_Delta"	 =>	"reductive",
		"Sat_NF"	 =>	"unspecific",
		"Sat_small"	 =>	"unspecific",
		"Sat_SOB"	 =>	"oxidative",
		"Sat_SOB2"	 =>	"oxidative",
		"Sat_SRO"	 =>	"reductive",
	); 	

}
######################################################################################################
sub Userinput{
	$in =$_[0];
	$fh = FileHandle->new;
	open($fh, "< $in") or die $!;
	while($line=<$fh>){
		chomp($line);
		@list=split(/\s+/, $line);
		$model=$list[0];
		$score=$list[1];
		$evalue=$list[2];
		$score{$model}=$score;
		$evalue{$model}=$evalue;
		$dir{$model}="not_defined";

	}
	close($fh);
	my @p=("AprA_arch", "AprA_ass", "AprA_Chlorobi", "AprA_SOB_M", "AprA_SOB_Q", "AprA_SRO", "AprB_arch", "AprB_ass", "AprB_Chlorobi", "AprB_SOB_M", "AprB_SOB_Q", "AprB_SRO", "AprM", "DsrA_arch", "DsrA_SDM", "DsrA_SOB", "DsrA_SRO", "DsrB_arch", "DsrB_SDM", "DsrB_SOB", "DsrB_SRO", "DsrC_arch", "DsrC_SDM", "DsrC_SOB", "DsrC_SRO", "DsrD1_SDM", "DsrD1", "DsrD2", "DsrD3", "DsrD4", "DsrE", "DsrF", "DsrH", "DsrJ_Chlorobi", "DsrJ_SDM", "DsrJ_SOB", "DsrJ_SRO", "DsrK_Chlorobi", "DsrK_MK_SRO", "DsrK_SDM", "DsrK_SOB", "DsrK_SRO", "DsrL_SDM", "DsrL_SOB", "DsrM_Chlorobi", "DsrM_MK_SRO", "DsrM_SDM", "DsrM_SOB", "DsrM_SRO", "DsrO_Chlorobi", "DsrO_SDM", "DsrO_SOB", "DsrO_SRO", "DsrP_Chlorobi", "DsrP_SDM", "DsrP_SOB", "DsrP_SRO", "DsrT_Chlorobi", "DsrT_SRO", "HdrA_SOB", "HdrA_SRO", "HdrA1_methano", "HdrA2_methano", "HdrB_SOB", "HdrB_SRO", "HdrB1_methano", "HdrB2_methano", "HdrB2_SOB", "HdrC_SOB", "HdrC_SRO", "HdrC1_methano", "HdrC2_methano", "HdrC2_SOB", "QmoA_SOB", "QmoA_SRO", "QmoA2_SRO", "QmoB_SOB", "QmoB_SRO", "QmoB2_SRO", "QmoC_SOB", "QmoC_SRO", "Sat_arch", "Sat_ass", "Sat_CF", "Sat_Chlorobi", "Sat_Delta", "Sat_NF", "Sat_small", "Sat_SOB", "Sat_SOB2", "Sat_SRO");

	foreach $model (@p){
		if(!exists $score{$model}){
			die "$in\tis incomplete - threshold missing\n";
		}
	}

}


######################################################################################################
#
#		usage
#
######################################################################################################
sub usage {
my $program = basename($0);
print <<END_OF_USAGE;
DESCRIPTION:
   DiSCo: Dissimilatory Sulfur metabolism Classification tOol
	-> filtering step

USAGE:
\n./$program -f <filename> -s <sequence_infile> -t <outfile> \[-o <filename>\] \[-d\] \[-h\] \n\n
  -f HMMer_outfile
  -o outfile_sequence filename  
  -s sequence_infile
  -d s <number>\t\tseparator for hitlist
     0 = TAB=\\t\t\t[default: tab separated]
     1 = ,
     2 = ;
     3 = SPACE=" "
  -t seperated list of positive hits
  -c tab seperated table with user-defined cut-offs: column 1: model, column 2: score, column 3: E-value
  -h print usage\n\n

 example
  perl filter_DiSCo.pl -f example.DiSCo.txt -s example.faa -t example.DiSCO.filtered.txt
  perl filter_DiSCo.pl -f example.DiSCo.txt -s example.faa -t example.DiSCO.filtered.txt -s 1 -o example.DiSCO.filtered.faa

END_OF_USAGE
}
######################################################################################################
