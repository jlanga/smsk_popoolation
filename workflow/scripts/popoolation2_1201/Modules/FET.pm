{
	package FET;
	use strict;
	use warnings;
	use List::Util qw[min max];
	use strict;
        use warnings;
        use FindBin qw/$Bin/;
        use lib "$Bin";
	use Text::NSP::Measures::2D::Fisher::twotailed;
    
        require Exporter;
        our @ISA = qw(Exporter);
        our @EXPORT=qw(get_fetcalculator);
	
	
    
	sub get_fetcalculator{
	    
	    my $popcount=shift;
	    my $winSumMethod=shift;
	    
	    return sub
	    {
		my ($window)=@_;
		die "No window provided" unless $window;
            
		my $data=$window->{data};
		die "No data for window $window->{chr}" unless $data;
		my $snps= [grep {$_->{ispuresnp}} @$data];
	    
		FET::_calculatePvalue($snps,$popcount,$winSumMethod);
	    }
	}
	
	sub _calculatePvalue
	{
	    my ($data,$popcount,$winSumMethod) = @_;
    
	    my $toret = {};
	    my $snpcount=@$data; 
    
	    my $samp=[];
	    
	    #calculate the pairwise index
	    my $index = [];
	    for(my $i=0; $i<$popcount; $i++)
	    {
		for(my $k=$i+1;$k<$popcount; $k++)
		{
		    my $id=$i.":".$k;
		    push @$index,$id;
		}
	    }
    
	    my $allsnps = [];
	    foreach my $d (@$data) {
		    
		next unless $d->{ispuresnp};
		    
		my $refchar = $d->{refchar};	
		my $cosa=$d->{samples}; #countsamples
		
		my ($major,$minor) = MajorAlleles::get_major_minor_alleles($cosa);
		
		my $snp = {};
		for(my $i=0; $i<@$cosa; $i++) {
		    my $t=$cosa->[$i];
		    my $cov=$t->{eucov};
		    $snp->{$i}=$t;
		}
		
		my $snp_pvalue = FET::_get_snpwise_pvalue($index,$snp,$major,$minor,$popcount);
		push @$allsnps,$snp_pvalue;
	    }
	    
	    $snpcount=@$allsnps;
	    
	    if ($snpcount==0) {
		for my $id (@$index) {
		    
		    $toret->{$id}="na";
		    
		}
	    }
	    else {
		
		for my $id (@$index) {
		    
		    my @pvalues=();
		    foreach my $snp (@$allsnps)
		    {
			push @pvalues,$snp->{$id};
		    }
		    
		    my $product_pvalue;
		    if($winSumMethod eq "multiply")
		    {
			$product_pvalue = multiply(\@pvalues);
		    }
		    elsif($winSumMethod eq "median")
		    {
			$product_pvalue = median(\@pvalues);
		    }
		    elsif($winSumMethod eq "geometricmean")
		    {
			$product_pvalue = geometricmean(\@pvalues);
		    }
		    else
		    {
			die "unrecognised method $winSumMethod";
		    }
		    
		    $toret->{$id}= $product_pvalue;
		}
	    }
	    
	    return $toret;
	}
	
	
	sub multiply
	{
		my $pvalues=shift;
		return "na" unless @$pvalues;
		my $toret=0;
		foreach my $pv (@$pvalues)
		{
			if($pv>0)
			{
				my $tmp=log($pv)/log(10);	
				$toret+=$tmp;
			}
			else
			{
				return "inf";
			}
		}
		$toret*=-1 if $toret<0;
		return $toret;
	}
	
	sub geometricmean
	{
		my $pvalues=shift;
		return "na" unless @$pvalues;
		my $toret=0;
		foreach my $pv (@$pvalues)
		{
			if($pv>0)
			{
				my $tmp=log($pv)/log(10);	
				$toret+=$tmp;				
			}
			else
			{
				return "inf";
			}
		}
		my $count=@$pvalues;
		$toret=$toret/$count;
		$toret*=-1 if $toret<0;
		return $toret;
	}
	
	sub median
	{
		my $pvalues=shift;
		return "na" unless @$pvalues;
		
		$pvalues=[sort {$a<=>$b} @$pvalues];
		my $count=@$pvalues;
		my $index=0;
		my $midcount=int($count/2);
		if($count%2==0)
		{
			my @tp=();
			push @tp,$pvalues->[$midcount-1];
			push @tp,$pvalues->[$midcount];
			return geometricmean(\@tp);
		}
		else
		{
			my @tp=();
			push @tp,$pvalues->[$midcount];
			return geometricmean(\@tp);
		}
	}

## ori
	#    else {
	#	
	#	for my $id (@$index) {
	#	    
	#	    my $product_pvalue=0;
	#	    
	#	    foreach my $snp (@$allsnps) {
	#		#print "$id:$snp->{$id}\t";	
	#		last if $product_pvalue eq "inf";
	#		
	#		my $pval=$snp->{$id};
	#	
	#		if($pval>0)
	#		{
	#			my $temp=log($pval)/log(10);
	#			$product_pvalue+=$temp;
	#		}
	#		else
	#		{
	#			$product_pvalue="inf";
	#		}
	#	    }
	#	    
	#	    $product_pvalue = $product_pvalue * -1 if $product_pvalue <0;
	#	    $toret->{$id}= $product_pvalue;
	#	}
	#    }
	#    
	#    return $toret;
	#}
	
	
	
	sub _get_snpwise_pvalue {
	# for a given SNP, a given major allele and a given minor allele
	# a specified pairwise comparision
    
	    my ($index,$snp,$major,$minor,$popcount) = @_;
	    
	    my $snp_pvalue = {};
	    foreach my $id (@$index) {
		
                my ($major1,$major2,$minor1,$minor2) = (0,0,0,0);
		my $string="";
		my $pvalue=1;

		my ($a,$b) = split('\:',$id);
		
                ($major1,$major2,$minor1,$minor2) = ($snp->{$a}->{$major},$snp->{$b}->{$major},$snp->{$a}->{$minor},$snp->{$b}->{$minor});

                $pvalue = FET::_get_Fisher_twotailed_pvalue($major1,$major2,$minor1,$minor2);
		$snp_pvalue->{$id}=$pvalue;
	    }
	    
	    return $snp_pvalue;
	    
	}
        
        sub _get_Fisher_twotailed_pvalue {
            my ($major1,$major2,$minor1,$minor2) = @_;
            
 
            
            my $np1=0;
            my $np2=0;
            my $n1p=0;
            my $n2p=0;
            my $npp=0;
            
            $np1 = $major1+$major2;
            $np2 = $minor1+$minor2;
            $n1p = $major1+$minor1;
            $n2p = $major2+$minor2;
            $npp = $n1p+$n2p;
            
            my $twotailed_value = calculateStatistic( n11=>$major1,
                                          n1p=>$n1p,
                                          np1=>$np1,
                                          npp=>$npp);
            
            my $pvalue=1;
            my $errorCode;
            if( ($errorCode = getErrorCode()))
            {
              $pvalue=1;
            }
            else
            {
              $pvalue = $twotailed_value;
            }
            
            return $pvalue;
    
        
        }

	
}