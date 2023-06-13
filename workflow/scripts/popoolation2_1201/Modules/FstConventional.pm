{
    package FstConventional;
    use strict;
    use warnings;
    use FindBin qw/$Bin/;
    use lib "$Bin";
    use List::Util qw[min max];
    
    require Exporter;
    our @ISA = qw(Exporter);
    our @EXPORT=qw(get_conventional_fstcalculator);
    
    
    sub get_conventional_fstcalculator
    {
        my $poolar=shift;
        my $popcount=shift;
        
        return sub
        {
            my $window=shift;
            die "No window provided" unless $window;
            
            my $data=$window->{data};
            die "No data for window $window->{chr}" unless $data;
            my $snps= [grep {$_->{ispuresnp}} @$data];
            
            my $frequencies=_calculateSNPFrequencies($snps,$popcount);
            my $hetero=_calculatePivalues($frequencies,$poolar);
            my $fstvalues=_calculateFstValues($hetero);
            return $fstvalues;
        }
    }




    sub _calculateSNPFrequencies
    {
        my $data=shift;
        my $popcount=shift;
        #my $popcount=shift;
        my $snpcount=@$data; 

        my $samp=[];
        

        foreach my $d (@$data)
        {
            next unless $d->{ispuresnp};
            my $cosa=$d->{samples}; #countsamples
            
            die "sample sizes are not agreeing" unless @$cosa==$popcount;
            
            my $freqsa=[];
            for(my $i=0; $i<@$cosa; $i++)
            {
                my $t=$cosa->[$i];
                my $cov=$t->{eucov};
                my $e={
                    index=>$t->{index},
                    eucov=>$cov,
                    A=>$t->{A}/$cov,
                    T=>$t->{T}/$cov,
                    C=>$t->{C}/$cov,
                    G=>$t->{G}/$cov
                };
                $freqsa->[$i]=$e;
            }
            push @$samp,$freqsa;
        }
        
        
        # reorganize the strucure
        my $toret={};
        # first initiate the thing
        for(my $i=0; $i<$popcount; $i++)
        {
            $toret->{$i}=[];
        }
        
        for(my $i=0; $i<@$samp; $i++)
        {
            my $s=$samp->[$i];
            for(my $k=0; $k<$popcount; $k++)
            {
                push @{$toret->{$k}},$s->[$k];
            }
        }
        
        #calculate the pairwise frequencies
        for(my $i=0; $i<$popcount; $i++)
        {
            for(my $k=$i+1;$k<$popcount; $k++)
            {
                my $id=$i.":".$k;
                my $ar=[];
                for(my $count=0; $count<$snpcount; $count++)
                {
                    my $a=$samp->[$count][$i];
                    my $b=$samp->[$count][$k];
                    push @$ar,
                    {
                        eucov=>min($a->{eucov},$b->{eucov}),####################
                        A=>($a->{A}+$b->{A})/2,
                        T=>($a->{T}+$b->{T})/2,
                        C=>($a->{C}+$b->{C})/2,
                        G=>($a->{G}+$b->{G})/2
                    }
                }
                $toret->{$id}=$ar; 
            }
        }
        
        my $control;
        while(my($key,$value)=each(%$toret))
        {
            $control=@$value unless $control;
            die "number of snps in $key is not in agreement with the others" unless $control == @$value;
        }
        
        return ($toret);
    }
    
    sub _calculatePivalues
    {
        my $frequencies=shift;
        my $poolar=shift;
        
        my $toret={};
        
        while(my($key,$value)=each(%$frequencies))
        {
            my $poolsize=0;
            
            if($key=~/:/)
            {
                my($a,$b)=split/:/,$key;
    
                $poolsize=min($poolar->[$a],$poolar->[$b]);#####################
            }
            else
            {
                $poolsize=$poolar->[$key];
            }
            
            my $pi=0;
            
            $pi=_pi($value,$poolsize);
            $toret->{$key}=$pi;
        }
        return ($toret);
    }
    
    sub _pi
    {

        my $snps=shift; # a simple array of snps
        my $poolSize=shift;
    
        my $pi_sum=0;
        foreach my $snp (@$snps)
        {
            my $pi_snp = _uncorrectedPiPerSNPFromFreqs($snp);
            $pi_sum += $pi_snp;
        }
        $pi_sum*=($poolSize/($poolSize-1));
        return $pi_sum;
    }


    sub _uncorrectedPiPerSNPFromFreqs
    {
        my $snp = shift;
        
        my $M=$snp->{eucov};
        my $pi_snp=1;
        $pi_snp-=($snp->{A})**2;
        $pi_snp-=($snp->{T})**2;
        $pi_snp-=($snp->{C})**2;
        $pi_snp-=($snp->{G})**2;
        
        #this term is indepent of the correction
        $pi_snp*=($M/($M-1));
        return $pi_snp;
    }
    
    
    sub  _calculateFstValues
    {
        my $hetero=shift;
        
        my $fsthash={};
        while(my($key,$value)=each(%$hetero))
        {
            my $poolsize=0;
            
            if($key=~/:/)
            {
                my($a,$b)=split/:/,$key;
                
                my $h1=$hetero->{$a};
                my $h2=$hetero->{$b};
                my $av=($h1+$h2)/2;
                my $h12=$value;
                
                my $val="0";
                $val=(($h12-$av)/$h12) if $h12;
                
                $fsthash->{$key}=$val;
            }
        }
        return $fsthash;
    }

    

}
1;