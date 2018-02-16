use Scalar::Util qw(looks_like_number);
open OLD, "PriorMed_training.csv" or die;
while ($line=<OLD>){
	chomp $line;
	@t=split ',', $line;
	$id=$t[2];
	$i=14;
	while ($i<17){
		if (looks_like_number($t[$i])){
			if ($line=~/YES/){
				if ($t[$i+1]=~/DAYS/){
					if ($t[$i]>$ref{$id}){
						$ref{$id}=$t[$i];
					}
				}
			}
		}
			
	
		$i++;
	}
}
@all=keys %ref;
$t=0;
$c=0;
foreach $aaa (@all){
	$t+=$ref{$aaa};
	$c++;
}
$avg=$t/$c;

open OLD, "CoreTable_training.csv" or die;
while ($line=<OLD>){
	chomp $line;
	if ($line=~/VEN/){
		@table=split ",", $line;	
		$name=$table[2];
		if (exists $ref{$name}){
			print "$name\t$table[3]\t$table[4]\t$ref{$name}\n";
		}else{
			print "$name\t$table[3]\t$table[4]\t$avg\n";
		}
	}
}

