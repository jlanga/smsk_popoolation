#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Path;
use File::Basename; # to get the file path, file name and file extension
use FindBin qw/$RealBin/;
use lib "$RealBin/Modules";
use Test;
use Pileup;
use SyncSlider;
use FstConventional;
use FstAsymptUnbiased;
use MaxCoverage;


# --pool-size 500 --min-count 2 --min-coverage 4 --window-size 1000 --step-size 1000 --input test/snp.merge --output test/test.fst
# --pool-size 500 --min-count 2 --min-coverage 4 --window-size 1000 --step-size 1000 --input /Users/robertkofler/dev/PopGenTools/test/trmod.sync --output /Users/robertkofler/dev/PopGenTools/test/test.fst 

my $input;
my $output="";
my $help=0;
my $test=0;

my $windowsize=1000;
my $step=1000;

my $mincount=2;
my $mincoverage=4;
my $poolsize;
my $minCoverageFraction=0.0;
my $usermaxcoverage;
my $asunbiased=0;
my $suppressna=0;
#my $correction="VarianceExactCorrection";


GetOptions(
    "input=s"	        =>\$input,
    "output=s"	        =>\$output,
    "min-count=i"       =>\$mincount,
    "min-coverage=i"    =>\$mincoverage,
    "max-coverage=s"    =>\$usermaxcoverage,
    "min-covered-fraction=f"    =>\$minCoverageFraction,
    "pool-size=s"       =>\$poolsize,
    "window-size=i"     =>\$windowsize,
    "step-size=i"       =>\$step,
    "karlsson-fst"      =>\$asunbiased,
    "suppress-noninformative"       =>\$suppressna,
    "test"              =>\$test,
    "help"	        =>\$help
) or pod2usage(-msg=>"Wrong options",-verbose=>1);

# too many arguments should result in an error
pod2usage(-msg=>"Wrong options",-verbose=>1) if @ARGV;
pod2usage(-verbose=>2) if $help;
FstTest::runTests() if $test;
pod2usage(-msg=>"Input file does not exist",-verbose=>1) unless -e $input;
pod2usage(-msg=>"No output file has been provided",-verbose=>1) unless $output;
pod2usage(-msg=>"Minimum coverage <1 not allowed",-verbose=>1) if $mincoverage<1;
pod2usage(-msg=>"Poolsize has to be provided",-verbose=>1) unless $poolsize;
pod2usage(-msg=>"Minimum coverage must be equal or larger than minimum count",-verbose=>1) unless $mincoverage>= $mincount;
pod2usage(-msg=>"Maximum coverage has to be provided",-verbose=>1) unless $usermaxcoverage;

my $paramfile=$output.".params";
open my $pfh, ">",$paramfile or die "Could not open $paramfile\n";
print $pfh "Using input\t$input\n";
print $pfh "Using output\t$output\n";
print $pfh "Using min-count\t$mincount\n";
print $pfh "Using min-coverage\t$mincoverage\n";
print $pfh "Using max-coverage\t$usermaxcoverage\n";
print $pfh "Using min-covered-fraction\t$minCoverageFraction\n";
print $pfh "Using pool-size\t$poolsize\n";
print $pfh "Using window-size\t$windowsize\n";
print $pfh "Using step-size\t$step\n";
print $pfh "Using asympt-unbiased\t$asunbiased\n";
print $pfh "Using suppress na\t$suppressna\n";
print $pfh "Using test\t$test\n";
print $pfh "Using help\t$help\n";
close $pfh;

my $maxcoverage=get_max_coverage($input,$usermaxcoverage);
open my $ofh, ">$output" or die "Could not open output file";
my $reader=SyncSlider->new($input,$windowsize,$step,$mincount,$mincoverage,$maxcoverage);

# file count
my $popcount=$reader->count_samples();

# handle the poolsize
my $poolAr=Utility::_getPoolSize($poolsize,$popcount);
 
my $fstCalculator;
if($asunbiased)
{
    $fstCalculator=get_asymptunbiased_fstcalculator($popcount);
}
else
{
    $fstCalculator=get_conventional_fstcalculator($poolAr,$popcount);
}



while(my $window=$reader->nextWindow())
{
        my $chr=$window->{chr};
        my $pos=$window->{middle};
        my $win=$window->{window};
        my $above=$window->{count_covered};
        my $snpcount=$window->{countpuresnp};
        my $avmincov=$window->{avmincov};
        my $data=$window->{data};
        next unless @$data;
        my $coveredFrac=$above/$win;
        
        my $fsth=$fstCalculator->($window);
        
        my $sufficientCovered = $coveredFrac>=$minCoverageFraction;
        next if(not $sufficientCovered and $suppressna);
        next if(not $snpcount and $suppressna);
        
        $coveredFrac=sprintf("%.3f",$coveredFrac);
        $avmincov=sprintf("%.1f",$avmincov);
        
        my $str=Utility::formatOutput($fsth,$popcount,$sufficientCovered);

        print $ofh "$chr\t$pos\t$snpcount\t$coveredFrac\t$avmincov\t$str\n";
        

}
close $ofh;




exit;




{
    package Utility;
    use strict;
    use warnings;
    use List::Util qw[min max];
    use FindBin qw/$RealBin/;
    use lib "$RealBin/Modules";


    
    
    sub formatOutput
    {
        my $fsth=shift;
        my $popcount=shift;
        my $sufcovered=shift;
        
        
        my $e=[];
        for(my $i=0; $i<$popcount; $i++)
        {
            for(my $k=$i+1; $k<$popcount; $k++)
            {
                my $key=$i.":".$k;
                my $val=$fsth->{$key};
                $val="na" unless $sufcovered;
                
                $val=sprintf("%.8f",$val) unless $val eq "na";
                my $str=($i+1).":".($k+1)."=".$val;
                push @$e,$str;
            }
        }
        my $toret=join("\t",@$e);
        return $toret;
    }
    
    
 
    


    
   
    

    
    sub _getPoolSize
    {
        my $poolsize=shift;
        my $popCount=shift;
        my $poolAr=[];
        if($poolsize=~/:/)
        {
            $poolAr=[split /:/, $poolsize];
        }
        else
        {
            for(1..$popCount)
            {
                push @$poolAr,$poolsize;
            }
        }
        die "provide poolsize does not fit with input file" if scalar(@$poolAr) != $popCount;
        return $poolAr;
    }
 
}

{
    package FstTest;
    use FindBin qw/$RealBin/;
    use lib "$RealBin/Modules";
    use Test;
    use Pileup;
    use Test::TSyncSlider;
    use Test::TFstCalculation;
    use Test::TSynchronized;
    use Test::TMaxCoverage;
    

    sub runTests
    {
        run_MaxCoverageTests();
        run_SynchronizedTests();
        run_SyncSliderTests();
        run_FstCalculationTests();
        exit;
    }
    
}

    #"input=s"	    =>\$input,
    #"output=s"	    =>\$output,
    #"min-count=i"   =>\$mincount,
    #"min-coverage=i"=>\$mincoverage,
    #"pool-size=i"   =>\$poolsize,
    #"window-unit=s"  =>\$windowunit,
    #"window-size=i"  =>\$windowsize,
    #"step-size=i"   =>\$step,
    #"test"          =>\$test,
    #"help"	    =>\$help

=head1 NAME

fst-sliding.pl - Calculate pairwise Fst-values for a set of populations

=head1 SYNOPSIS

 fst-sliding.pl --input populations.merged --output population.fst --min-count 2 --min-coverage 4 --max-coverage 2% --pool-size 500 --window-size 1000 --step-size 1000

=head1 OPTIONS

=over 4

=item B<--input>

The input file. Has to be synchronized pileup file. Mandatory parameter

=item B<--output>

The output file. Mandatory parameter

=item B<--karlsson-fst>

flag; switch to an alternative method for calculating the pairwise fst which was introduced by Karlsson et al. (2007), use C<--help> for more information; default=not set

=item B<--min-count>

the minimum count of the minor allele. used for SNP identification. SNPs will be identified considering all populations simultanously. default=2

=item B<--min-coverage>

the minimum coverage; used for SNP identification, the coverage in ALL populations has to be higher or equal to this threshold, otherwise no SNP will be called. default=4

=item B<--max-coverage>

The maximum coverage; All populations are required to have coverages lower or equal than the maximum coverage; Mandatory
The maximum coverage may be provided as one of the following:

 '500' a maximum coverage of 500 will be used for all populations
 '300,400,500' a maximum coverage of 300 will be used for the first population, a maximum coverage of 400 for the second population and so on
 '2%' the 2% highest coverages will be ignored, this value is independently estimated for every population

=item B<--min-covered-fraction>

the minimum fraction of a window being between min-coverage and max-coverage in ALL populations; float; default=0.0

=item B<--window-size>

the size of the sliding window. Measured in C<--window-unit>; default=1000

=item B<--step-size>

the size of the sliding window steps. Measured in C<--window-unit>; default=1000

=item B<--pool-size>

the size of the population pools; May be provided for each population individually; mandatory parameter

 --pool-size 500 .. all populations have a pool size of 500
 --pool-size 500:300:600 .. first population has a pool size of 500, the seccond of 300 etc;
   the number of pools has to fit with the number of populations provided in the file

=item B<--suppress-noninformative>

Suppress output for windows with no SNPs or insufficient coverage; default=off

=item B<--test>

Run the unit tests for this script. 

=item B<--help>

Display help for this script

=back

=head1 DETAILS

Combine different pileup files

=head2 INPUT

A synchronized file, for example

 Unknown_group_104	5943	N	0:0:10:0:0:0	0:0:10:0:0:0	0:0:10:0:0:0
 Unknown_group_104	5944	N	0:8:0:0:0:0	0:8:0:0:0:0	0:8:0:0:0:0
 Unknown_group_104	5945	N	0:0:0:8:0:0	0:0:0:8:0:0	0:0:0:8:0:0
 Unknown_group_104	5946	N	0:0:9:0:0:0	0:0:9:0:0:0	0:0:9:0:0:0
 Unknown_group_104	5947	N	0:7:0:0:0:0	0:7:0:0:0:0	0:7:0:0:0:0
 
 col1: reference contig (chromosome)
 col2: position in the reference contig
 col3: reference character
 col4: population 1
 col5: population 2
 coln: population n
 
 population data are in the form
 A:T:C:G:N:*
 A: count of character A
 T: count of character T
 C: count of character C
 G: count of character G
 N: count of character N
 *: deletion, count of deletion

Only the characters A,T,C,G are considered for the coverage.
The SNP site is ignored if a deletion is found in any population

=head2 OUTPUT

The output will be given for every sliding window, e.g.:

 2L      68500   360     1.000   62.1    1:2=0.01873725  1:3=0.02131245  2:3=0.01521177
 2L      69500   118     1.000   71.9    1:2=0.00969479  1:3=0.00116059  2:3=0.00905794
 2L      70500   269     1.000   63.6    1:2=0.01955417  1:3=0.01547995  2:3=0.01300569

 col1: reference contig (chromosome)
 col2: mean position of the sliding window
 col3: number of SNPs found in the window (not considering sites with a deletion) 
 col4: fraction of the window which has a sufficient coverage (min. coverage <= cov <= max. coverage) in every population;
 col5: average minimum coverage in all populations

 ....
 1:2 the pairwise Fst for population 1 and 2
 1:3 the pairwise Fst for population 1 and 3

Note to col5: first the minimum coverage for every position is recorded in all populations;
This min coverage is then subsequently averaged over the whole windows.
It may, for example, happen that the min coverage of position 2L:5 is from population 1 but the min coverage of position 2L:6 from population 2.
Only positions which are sufficiently covered are considered.

=head2 Fst calculation

By default Fst is calculated from the allele-frequencies (not from the allele-counts) using the standard equation as shown in "Hartl and Clark (2007): Principles of Population Genetics"
 
 Fst = (Pi_total - Pi_within) / Pi_total
 Pi_within = (Pi_population1 + Pi_population2)/ 2
 Pi: 1 - fA ^ 2 - fT ^ 2 -fC ^ 2 - fG ^ 2
 fN:  frequency of nucleotide N
 Pi_total: for the total Pi the allele frequencies of the two
     populations are averaged and Pi is calculated as shown above

if the alternative method C<--karlsson-fst> is used, Fst is calculated  from the allele-counts according to:
Karlsson et al. (2007): Efficient mapping of mendelian traits in dogs through genome-wide association, Nature genetics 


=head1 AUTHORS

Robert Kofler

Christian Schloetterer

=cut

