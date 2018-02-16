		$j=2;
		while ($j<10){
				print "$j\t";
				system "perl ~/source/data_stat/cor_column.pl max.ref 2 test_all.dat $j";
			$j++;
		}
		
