#!/usr/bin/perl

use Scalar::Util qw(looks_like_number);
open FILE, "PriorMed_leaderboard.csv" or die;
while ($line=<FILE>){
	if ($line=~/HORMONAL THERAPY/){
		chomp $line;
		@table=split ',', $line;
		$n=scalar (@table);
		$i=36;
		while ($i<$n){
			if (looks_like_number($table[$i])){
				if ((looks_like_number($table[$i+1]))||($table[$i+1] eq '.')){
					
					if ($table[$i+1] eq '.'){

						$table[$i+1]=0;
					}
					$start=$table[$i];
					$end=$table[$i+1];
					if (exists $record{$table[2]}){

						$record{$table[2]}.="\t$start,$end";
					}else{
						$record{$table[2]}="$start,$end";
					}
				}
			}
				
			$i++;
		}
	}
}
close FILE;


@all=keys %record;

foreach $aaa (@all){
	@list=split "\t", $record{$aaa};
	$min=0;
	$max=0;
	foreach $lll (@list){
		@t=split ',', $lll;
		$val=$t[1]-$t[0];
		if ($val>$ref{$aaa}){
			$ref{$aaa}=$val;
		}
	}
}


@all=keys %ref;
@total=();
$t=0;
$c=0;
foreach $aaa (@all){
	push @total, $ref{$aaa};
	$t+=$ref{$aaa};
	$c++;
}
@total=sort{$a<=>$b}@total;
$n=int(scalar(@total)/2);
$avg=$total[$n];
$t=$t/$c;




open OLD, "CoreTable_leaderboard_new.csv" or die;
open NEW, ">max.ref" or die;
<OLD>;
while ($line=<OLD>){
	chomp $line;
	@table=split ",", $line;
	print NEW "$table[2]\t";
	if (exists $ref{$table[2]}){
		print NEW "$ref{$table[2]}";
	}else{
		print NEW "$t";
		print "$table[2]\t$t\n";
	}
	print NEW "\n";
}
close OLD;
close NEW;

	
		
	
