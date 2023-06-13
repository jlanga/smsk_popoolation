#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use FindBin qw/$RealBin/;
use lib "$RealBin/Modules";
use FET;
use MaxCoverage;
use SyncSlider;
use SynchronizeUtility;
use MajorAlleles; # get the two major allele
use Test;
use List::Util qw[min max];

# Author: Ram Vinay Pandey 

# Define the variables
my $input;
my $output="";
my $help=0;
my $test=0;
my $verbose=1;

my $winSumMethod="multiply";
my $windowsize=1;
my $step=1;
my $mincount=2;
my $mincoverage=4;
my $usermaxcoverage;
my $minCoverageFraction=0.0;
my $suppressna=0;

# -input /Volumes/Main/popoolation2/test.syn -output /Volumes/Main/popoolation2/test.out -min-count 2 -min-coverage 10 -max-coverage 200,300,200,500 -window-size 1 -step-size 1 -min-covered-fraction 0.6    
    
GetOptions(
    "input=s"	    			=>\$input,
    "output=s"	    			=>\$output,
    "min-count=i"   			=>\$mincount,
    "min-coverage=i"			=>\$mincoverage,
    "max-coverage=s"			=>\$usermaxcoverage,
    "window-size=i"  			=>\$windowsize,
    "step-size=i"   			=>\$step,
    "min-covered-fraction=f"		=>\$minCoverageFraction,
    "window-summary-method=s"		=>\$winSumMethod,
    "suppress-noninformative"       	=>\$suppressna,
    "test"          			=>\$test,
    "help"	    			=>\$help
) or pod2usage(-msg=>"Wrong options",-verbose=>1);

# too many arguments should result in an error
pod2usage(-msg=>"Wrong options",-verbose=>1) if @ARGV;
pod2usage(-verbose=>2) if $help;
FETTest::runTests() && exit if $test;

$winSumMethod=lc($winSumMethod);
pod2usage(-msg=>"Invalid Window Summary method") if($winSumMethod ne "multiply" and $winSumMethod ne "geometricmean" and $winSumMethod ne "median");
pod2usage(-msg=>"Input file does not exist",-verbose=>1) unless -e $input;
pod2usage(-msg=>"No output file has been provided",-verbose=>1) unless $output;
pod2usage(-msg=>"Minimum coverage <1 not allowed",-verbose=>1) if $mincoverage<1;
pod2usage(-msg=>"Minimum coverage must be equal or larger than minimum count",-verbose=>1) unless $mincoverage>= $mincount;
pod2usage(-msg=>"Maximum coverage has to be provided",-verbose=>1) unless $usermaxcoverage;

my $paramfile=$output.".params";
open my $pfh, ">",$paramfile or die "Could not open $paramfile\n";
print $pfh "Using input\t$input\n";
print $pfh "Using output\t$output\n";
print $pfh "Using min-count\t$mincount\n";
print $pfh "Using min-coverage\t$mincoverage\n";
print $pfh "Using max-coverage\t$usermaxcoverage\n";
print $pfh "Using window-size\t$windowsize\n";
print $pfh "Using step-size\t$step\n";
print $pfh "using window summary method\t$winSumMethod\n";
print $pfh "Using min-covered-fraction\t$minCoverageFraction\n";
print $pfh "Using test\t$test\n";
print $pfh "Using help\t$help\n";
close $pfh;



open my $ofh, ">$output" or die "Could not open output file";

my $maxcoverage=get_max_coverage($input,$usermaxcoverage);
my $reader;
$reader=SyncSlider->new($input,$windowsize,$step,$mincount,$mincoverage,$maxcoverage);

# file count
my $popcount=$reader->count_samples();

my $fetCalculator;
$fetCalculator=FET::get_fetcalculator($popcount,$winSumMethod);


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
        
        my $sufficientCovered=$coveredFrac>=$minCoverageFraction;
        next if(not $sufficientCovered and $suppressna);
        next if(not $snpcount and $suppressna);
        
	my $feth=$fetCalculator->($window);
	
        $coveredFrac=sprintf("%.3f",$coveredFrac);
        $avmincov=sprintf("%.1f",$avmincov);
	
	my $str=Utility::formatOutput($feth,$popcount,$sufficientCovered);

        print $ofh "$chr\t$pos\t$snpcount\t$coveredFrac\t$avmincov\t$str\n";
	#print "$chr\t$pos\t$snpcount\t$coveredFrac\t$avmincov\t$str\n";
	

}


exit(0);





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
                
		unless($val eq "na" or $val eq "inf")
		{
		    $val=sprintf("%.8f",$val);    
		}
                
                my $str=($i+1).":".($k+1)."=".$val;
                push @$e,$str;
            }
        }
        my $toret=join("\t",@$e);
        return $toret;
    }
    
    
}


{
    package FETTest;
    use FindBin qw/$RealBin/;
    use lib "$RealBin/Modules";
    use Test;
    use Pileup;
    use Test::TSyncSlider;
    use Test::TFisher;
    use Test::TSynchronized;
    use Test::TMaxCoverage;
    

    sub runTests
    {
        run_MaxCoverageTests();
        run_SynchronizedTests();
        run_SyncSliderTests();
        run_FisherTests();
        exit;
    }
    
}


    #"input=s"	    			=>\$input,
    #"output=s"	    			=>\$output,
    #"min-count=i"   			=>\$mincount,
    #"min-coverage=i"			=>\$mincoverage,
    #"max-coverage=s"			=>\$usermaxcoverage,
    #"window-unit=s"  			=>\$windowunit,
    #"window-size=i"  			=>\$windowsize,
    #"step-size=i"   			=>\$step,
    #"min-covered-fraction=f"		=>\$minCoverageFraction,
    #"allele=i"   			=>\$allele,
    #"suppress-noninformative"       	=>\$suppressna,
    #"test"          			=>\$test,
    #"help"	    			=>\$help
    #
    
=head1 NAME

fisher-test.pl - Calculate pairwise p-values for a set of populations 

=head1 SYNOPSIS

 perl fisher-test.pl --input input.sync --output output.fet --min-count 2 --min-coverage 4 --max-coverage 2% 

=head1 OPTIONS

=over 4

=item B<--input>

The input file. Has to be synchronized file. Mandatory parameter

=item B<--output>

The output file. Mandatory parameter

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

the size of the sliding window. Measured in C<--window-unit>; default=1

=item B<--step-size>

the size of the sliding window steps. Measured in C<--window-unit>; default=1

=item B<--window-summary-method>

Specify the method by which the p-values of individual SNPs should be summarized; possible: geometricmean | median | multiply; default=multiply

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

 2L	5943	N	0:0:10:0:0:0	0:0:10:0:0:0	0:0:10:0:0:0
 2L	5944	N	0:8:0:0:0:0	0:8:0:0:0:0	0:8:0:0:0:0
 2L	5945	N	0:0:0:8:0:0	0:0:0:8:0:0	0:0:0:8:0:0
 2L	5946	N	0:0:9:0:0:0	0:0:9:0:0:0	0:0:9:0:0:0
 2L	5947	N	0:7:0:0:0:0	0:7:0:0:0:0	0:7:0:0:0:0
 
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

 4	60500	8	0.700	37.3	1:2=-6.10183008	1:3=-3.83272455	1:4=-6.62710556	2:3=-2.05551953	2:4=-3.12067424	3:4=-1.74211978
 4	61500	6	0.999	42.2	1:2=-2.02043203	1:3=-0.26377650	1:4=-4.45791884	2:3=-1.31524521	2:4=-2.44527364	3:4=-2.42537609
 4	62500	0	0.010	34.8	1:2=na	1:3=na	1:4=na	2:3=na	2:4=na	3:4=na
 4	63500	1	0.770	36.9	1:2=0.00000000	1:3=0.00000000	1:4=0.00000000	2:3=0.00000000	2:4=0.00000000	3:4=0.00000000
 4	64500	9	0.892	27.9	1:2=-6.49500602	1:3=-9.13530619	1:4=-3.76371334	2:3=-6.46589826	2:4=-9.22209785	3:4=-7.85951677

 col1: reference contig (chromosome)
 col2: mean position of the sliding window
 col3: number of SNPs found in the window (not considering sites with a deletion) 
 col4: fraction of the window which has a sufficient coverage (min. coverage <= cov <= max. coverage) in every population;
 col5: average minimum coverage in all populations

 ....
 1:2 -log(product of p-value of all SNPs in window) for comparision of population 1 and 2
 1:3 -log(product of p-value of all SNPs in window) for comparision of population 1 and 3
 1:4 -log(product of p-value of all SNPs in window) for comparision of population 1 and 4


=head1 AUTHORS

 Ram Vinay Pandey
 Robert Kofler
 Christian Schloetterer

=cut
