{
package MaxCoverage;
use strict;
use warnings;
use FindBin qw/$Bin/;
use lib "$Bin";
use Synchronized;
use SynchronizeUtility;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT     =qw(get_max_coverage);



sub get_max_coverage
    {

        my $syncfile=shift;
        my $usermaxcov=shift;
        
        my $popcount=get_popcount_forsyncfile($syncfile);
        
        my $maxcoverages=[];
        if($usermaxcov=~/%$/)
        {
            $usermaxcov=~s/%$//;
            $maxcoverages=_get_maxcov_bypercentoutlier($syncfile,$usermaxcov,$popcount);
        }
        elsif($usermaxcov=~/,/)
        {
            $maxcoverages=[split /,/ , $usermaxcov];
        }
        else
        {
            my @tempar=();
            for my $i (1..$popcount)
            {
                push @tempar,$usermaxcov;
            }
            $maxcoverages=\@tempar;
        }
        
        my $tempcount=@$maxcoverages;
        die "Number of populations $popcount does not agree with number of maximum coverages $tempcount" if($tempcount!=$popcount);
        
        return $maxcoverages;  
    }
    
    
    sub _get_maxcov_bypercentoutlier
    {
        my $syncfile=shift;
        my $pcoutlier=shift;
        my $popcount=shift;
        my $pp=get_basic_syncparser();
        my $minfraction=(100-$pcoutlier)/100;
        
        
        print "Computing maximum coverage cutoffs from the empirical distributions of coverages\n";
        open my $ifh,"<", $syncfile or die "Could not open input file";
        
        
        my $counter=[];
        # parse synchronized file and get the coverage distribution
        while(my $line=<$ifh>)
        {
            chomp $line;
            next unless $line;
            my $parsed=$pp->($line);
            my $samples=$parsed->{samples};
           
            for(my $i=0; $i<$popcount; $i++)
            {
                my $samp=$samples->[$i];
                my $cov=$samp->{eucov};
                next unless $cov; # ignore zero
                $counter->[$i][$cov]++;
            }
        }
        
        # determine the maximum coverages
        my $maxcoverages=[];
        for(my $i=0; $i<$popcount; $i++)
        {
            my $covs=$counter->[$i]||[];
            my $sum=0;
            
            for(my $k=0;$k<@$covs; $k++)
            {
                my $t=$covs->[$k];
                next unless $t;
                $sum+=$t;
            }
            
            my $curcounter=0;
            my $lastcoverage=0;
            COVERAGE: for(my $k=0; $k<@$covs; $k++)
            {

                my $t=$covs->[$k];
                next unless $t;
                $curcounter+=$t;
                my $fraction=$curcounter/$sum;
                if($fraction > $minfraction)
                {
                    $maxcoverages->[$i]=$lastcoverage;
                    last COVERAGE;
                }
                $lastcoverage=$k;
            }
            $maxcoverages->[$i]=0 unless $maxcoverages->[$i];
            
        }
        print "Finished identifying maximum coverages\n";
        my $maxstr=join(",",@$maxcoverages);
        print "Result: '--max-coverage $pcoutlier%' is equivalent to '--max-coverage $maxstr'\n";

        return $maxcoverages;
    }

}

1;