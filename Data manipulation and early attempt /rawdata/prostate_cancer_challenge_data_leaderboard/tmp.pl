$i=2;
while ($i<10){
		$j=2;
		while ($j<10){
			if ($i==$j){}else{
				print "$i\t$j\t";
				system "perl ~/source/data_stat/cor_column.pl test_all.dat $i test_all.dat $j";
			}
			$j++;
		}
		
	$i++
}

