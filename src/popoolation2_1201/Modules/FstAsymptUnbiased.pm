{
    package FstAsymptUnbiased;
    use strict;
    use warnings;
    use FindBin qw/$Bin/;
    use lib "$Bin";
    
    require Exporter;
    our @ISA = qw(Exporter);
    our @EXPORT=qw(get_asymptunbiased_fstcalculator);
    our @EXPORT_OK=qw(get_asymptunbiased_fstcalculator reorganizeSNPs calculate_nk_dk);
    
    sub get_asymptunbiased_fstcalculator
    {
        my $popcount=shift;
        
        return sub
        {
            my $window=shift;
            die "No window provided" unless $window;
            my $data=$window->{data};
            die "No data for window $window->{chr}" unless $data;
            my $snps= [grep {$_->{ispuresnp}} @$data];
            
            my $refsnp=reorganizeSNPs($snps,$popcount);
            
            my $fsts=getfsts($refsnp,$popcount);
        }
    }
    
    
    sub getfsts
    {
        my $snps=shift;
        my $popcount=shift;
        
        my $fsth={};
        
        for(my $i=0; $i<$popcount; $i++)
        {
            for(my $k=$i+1;$k<$popcount; $k++)
            {
                my $key="$i:$k";
                
                my $sumnk=0;
                my $sumdk=0;
                
                foreach my $sh (@$snps)
                {
                    my $snp=$sh->{$key};
                    
                    my($nk,$dk)=calculate_nk_dk($snp);
                    $sumnk+=$nk;
                    $sumdk+=$dk;
                }
                my $val="na";
                $val=$sumnk/$sumdk if $sumdk;
                $fsth->{$key}=$val;
            }
        }
        return $fsth; 
    }
    
    sub calculate_nk_dk
    {
        my $snp=shift;
        if($snp->{n1}<=1 or $snp->{n2}<=1)
        {
            return (0,0);
        }
        
        my $h1=calculate_h1($snp);
        my $h2=calculate_h2($snp);
                    
        my $nksquare=(($snp->{a1}/$snp->{n1})-($snp->{a2}/$snp->{n2}))**2;
        my $sub=($h1/$snp->{n1}+$h2/$snp->{n2});
        my $nk=$nksquare-$sub;
                    
        my  $dk=$nk+$h1+$h2;
        return ($nk,$dk);
    }
    
    sub calculate_h1
    {
        my $snp=shift;
        
        my $val=($snp->{a1}*($snp->{n1}-$snp->{a1}))/($snp->{n1}*($snp->{n1}-1))
        
    }
    
    sub calculate_h2
    {
        my $snp=shift;
        
        my $val=($snp->{a2}*($snp->{n2}-$snp->{a2}))/($snp->{n2}*($snp->{n2}-1))
        
    }
    
    
    sub reorganizeSNPs
    {
        my $snps=shift;
        my $popcount=shift;
        
        my $snpref=[];
        
        foreach my $snp (@$snps)
        {
            my $sh={};
            for(my $i=0; $i<$popcount; $i++)
            {
                for(my $k=$i+1; $k<$popcount; $k++)
                {
                    my $a=$snp->{samples}[$i];
                    my $b=$snp->{samples}[$k];
                    my $e=get_a1a2n1n2($a,$b);
                    my $key="$i:$k";
                    $sh->{$key}=$e;
                }
            }
            push @$snpref,$sh;
        }
        return $snpref;
    }
    
    
    sub get_a1a2n1n2
    {
        my $asample=shift;
        my $bsample=shift;
        
        my @nuc=qw/A T C G/;
        my @freq;
        for my $n (@nuc)
        {
            push @freq,
            {
                nuc=>$n,
                avfreq=>(($asample->{$n}/$asample->{eucov}) + ($bsample->{$n}/$bsample->{eucov}))/2
            };
        }
        
        @freq=sort {$b->{avfreq}<=>$a->{avfreq}} @freq;
        
        my $major=$freq[0]->{nuc};
        my $minor=$freq[1]->{nuc};
        
        my $a1=$asample->{$major};
        my $n1=$asample->{$major}+$asample->{$minor};
        my $a2=$bsample->{$major};
        my $n2=$bsample->{$major}+$bsample->{$minor};
        
        return {
            a1=>$a1,
            n1=>$n1,
            a2=>$a2,
            n2=>$n2
        };
    }
}
1;