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
				$val=-$table[$i];
				if ($val>$ref{$table[2]}){
					$ref{$table[2]}=$val;
				}
			}
				
			$i++;
		}
	}
}
close FILE;



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
open NEW, ">long.ref" or die;
<OLD>;
while ($line=<OLD>){
	chomp $line;
	@table=split ",", $line;
	print NEW "$table[2]\t";
	if (exists $ref{$table[2]}){
		print NEW "$ref{$table[2]}";
	}else{
		print NEW "$avg";
		print "$avg\n";
	}
	print NEW "\n";
}
close OLD;
close NEW;

	
		
	
