# !/usr/bin/perl
#use strict;
use Getopt::Long;
my %opt=();
GetOptions(\%opt,"p=s","l=s","s:s","e:s");
my $usage=<<USAGE;

Function: extract the primer sequences from the results of PRIMER3

Contact: Zhonghua Zhang, zhangzhonghua.caas@gmail.com

Usage: $0 -p primer3_result -l primer_list -s primer_seq -e for_e_PCR

USAGE
die $usage if (!defined $opt{p} ||!defined $opt{l});
open(PRIMER,$opt{p})||die "Cant open the primer file.\n";
open(LIST,">".$opt{l});
if(defined $opt{s})
{
	open(SEQ,">".$opt{s});
}
if(defined $opt{e})
{
	open(PCR,">".$opt{e});
}
$/="\n=\n";

while(<PRIMER>)
{
	my $name="";
	my $primer_l="";
	my $primer_r="";
	my $left_TM=0;
	my $right_TM=0;
	my $size=0;
	my @tmp=();
	/PRIMER_SEQUENCE_ID=(.+)/;
	$name=$1;
	@tmp=split(/\|/,$name);
	if(/PRIMER_LEFT_0_SEQUENCE=(\w+)/)
	{
		$primer_l=$1;
		/PRIMER_RIGHT_0_SEQUENCE=(\w+)/;
		$primer_r=$1;
		if(defined $opt{s}){
			print SEQ ">L_".$name."\n".$primer_l."\n";
			print SEQ ">R_".$name."\n".$primer_r."\n";
		}
		if(defined $opt{e}){
			print PCR $tmp[0],"\t",$primer_l,"\t",$primer_r,"\t","50-1500","\n";
		}
		/PRIMER_LEFT_0_TM=(\S+)/;
		$left_TM=$1;
		/PRIMER_RIGHT_0_TM=(\S+)/;
		$right_TM=$1;
		/PRIMER_PRODUCT_0_SIZE=(\d+)/;
		$size=$1;
		local $i=0;
		for($i=0;$i<@tmp;$i++)
		{
			print LIST $tmp[$i],"\t";
		}
		print LIST $left_TM,"\t",$right_TM,"\t",$size,"\t",$primer_l,"\t",$primer_r,"\n";
		#print OF $name."\t".$primer_l."\t",$primer_r,"\n";
	}
	else
	{
		next;
	}
	
}
close(PRIMER);close(LIST);close(SEQ);close(PCR);
