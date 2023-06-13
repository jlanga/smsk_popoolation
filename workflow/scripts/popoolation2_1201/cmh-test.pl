#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Basename;
use File::Path;
use FindBin qw/$RealBin/;
use lib "$RealBin/Modules";
use MaxCoverage;
use Synchronized;
use SynchronizeUtility;
use MajorAlleles; # get the two major allele
use Test;


# Author: Ram Vinay Pandey
# Author: Robert Kofler


# Define the variables
my $input;
my $output="";
my $userpopulation;
my $help=0;
my $test=0;
my $verbose=1;

my $mincount=2;
my $mincoverage=4;
my $usermaxcoverage;
my $minpvalue=1.0;
my $removetemp=0;

# --input /Users/robertkofler/pub/PoPoolation2/Walkthrough/demo-data/cmh/small-test.sync --output /Users/robertkofler/pub/PoPoolation2/Walkthrough/demo-data/cmh/small-test.cmh --population 1,2,3,4 --min-count 2 --min-coverage 4 --max-coverage 200

GetOptions(
    "input=s"	    =>\$input,
    "output=s"	    =>\$output,
    "min-count=i"   =>\$mincount,
    "min-coverage=i"=>\$mincoverage,
    "max-coverage=s"=>\$usermaxcoverage,
    "population=s"  =>\$userpopulation,
    "min-pvalue=f"  =>\$minpvalue,
    "remove-temp"   =>\$removetemp,
    "test"          =>\$test,
    "help"	    =>\$help
) or pod2usage(-msg=>"Wrong options",-verbose=>1);

pod2usage(-verbose=>2) if $help;
CMHTest::runTests() && exit if $test;
pod2usage(-msg=>"A input file has to be provided\n",-verbose=>1) unless -e $input;
pod2usage(-msg=>"A output file has to be provided\n",-verbose=>1) unless $output;
pod2usage(-msg=>"Minimum coverage must be equal or larger than minimum count",-verbose=>1) unless $mincoverage>= $mincount;
pod2usage(-msg=>"Maximum coverage has to be provided",-verbose=>1) unless $usermaxcoverage;
pod2usage(-msg=>"The pairwise comparisions have to be provided (--population)",-verbose=>1) unless $userpopulation;



################# write param file

my $paramfile=$output.".params";
open my $pfh, ">",$paramfile or die "Could not open $paramfile\n";
print $pfh "Using input\t$input\n";
print $pfh "Using output\t$output\n";
print $pfh "Using min-count\t$mincount\n";
print $pfh "Using min-coverage\t$mincoverage\n";
print $pfh "Using max-coverage\t$usermaxcoverage\n";
print $pfh "Using population\t$userpopulation\n";
print $pfh "Using min-pvalue\t$minpvalue\n";
print $pfh "Remove temporary files\t$removetemp\n";
print $pfh "Using test\t$test\n";
print $pfh "Using help\t$help\n";
close $pfh;

my $maxcoverage=get_max_coverage($input,$usermaxcoverage);
my $populations=CMHUtil::resolve_population($userpopulation);
my $syncparser=get_sumsnp_synparser($mincount,$mincoverage,$maxcoverage);

my $rinput=$output.".rin";
my $routput=$output.".rout";

print "Reading sync file and writing temporary R output file\n";
CMHUtil::write_Rinput($input,$rinput,$syncparser,$populations);

print "Calling R, to calculate the Cochran-Mantel-Haenszel test statistic\n";
system("R --vanilla --slave <$rinput >$routput");

print "Parsing R-output and writing output file\n";
CMHUtil::write_output($routput,$output,$minpvalue);

if($removetemp)
{
	print "Removing temporary files\n";
	unlink($rinput);
	unlink($routput);
}
print "Done\n";

exit(0);



{
	package CMHUtil;
	use strict;
	use warnings;
	use List::Util qw[min max];
	use FindBin qw/$RealBin/;
	use lib "$RealBin/Modules";
	use MaxCoverage;
	use Synchronized;
	use SynchronizeUtility;
	use MajorAlleles; # get the two major allele
	
	sub write_output
	{
		my $routput=shift;
		my $output=shift;
		my $minpvalue=shift;
		
		open my $ifh,"<", $routput or die "Could not open input file\n";
		open my $ofh,">",$output or die "Could not open output file\n";
		
		while(1)
		{
			#[1] "2R\t2296\tN\t90:10:0:0:0:0\t100:0:0:0:0:0\t100:0:0:0:0:0\t100:0:0:0:0:0"
			#[1] 0.003583457
			my $line=<$ifh>;
			last unless $line;
			my $pvalue=<$ifh>;
			chomp $line; chomp $pvalue;
			$line=~s/^\S+\s//;
			$line=~s/^"//;
			$line=~s/"$//;
			$line=~s/\\t/\t/g;
			$pvalue=~s/^\S+\s//;
			next if $pvalue> $minpvalue;
			print $ofh $line."\t".$pvalue."\n";
		}
		close $ofh;
		close $ifh;
	}
	
	sub resolve_population
	{
		my $userpopulation=shift;
		die "At least two pairwise comparisions need to be specified, e.g.: 1-2,3-4" unless $userpopulation=~m/,/;
		die "At least two pairwise comparisions need to be specified, e.g.: 1-2,3-4" unless $userpopulation=~m/-/;
		
		
		my $populations=[];
		my @temp=split /,/,$userpopulation;
		foreach my $t (@temp)
		{
			my @a=split /-/,$t;
			die "At least two pairwise comparisions need to be specified" if scalar(@a) !=2;
			push @$populations,$a[0];
			push @$populations,$a[1];
		}
		
		die "Pairwise comparisions must be an even number (user provided $userpopulation)" if(scalar(@$populations) % 2);
		return $populations;
	}
	
	sub write_Rinput
	{
		my $syncfile=shift;
		my $rinput=shift;
		my $syncparser=shift;
		my $populations=shift;
		
		my $third_dim =int(scalar(@$populations)/2);
		my $dim_str="c(2,2,$third_dim)";
		
		open my $ifh, "<", $syncfile or die "Could not open input file";
		open my $ofh, ">", $rinput or die "Could not open routput file";
		while(my $line=<$ifh>)
		{
			chomp $line;
			my $e=$syncparser->($line);
			next unless $e->{ispuresnp};
			
			my $pop_str=_get_populationstring($e->{samples},$populations);
			my $ar_str="array($pop_str,dim=$dim_str)";
			my $mantel_str="mantelhaen.test($ar_str,alternative=c(\"two.sided\"))\$p.value";
			print $ofh "print(\"$line\")\n";
			print $ofh $mantel_str."\n";
		}
		close $ofh;
		close $ifh;
	}
	
	
	sub _get_populationstring
	{
		my $samples=shift;
		my $populations=shift;
		
		my ($major,$minor) = MajorAlleles::get_major_minor_alleles($samples);
		my @ar=();
		for(my $i=1; $i<@$populations; $i+=2)
		{
			my $basenr=$populations->[$i-1];
			my $derivednr=$populations->[$i];
			$basenr--; $derivednr--;
			my $base=$samples->[$basenr];
			my $derived=$samples->[$derivednr];
			
			push @ar,$base->{$major};
			push @ar,$derived->{$major};
			push @ar,$base->{$minor};
			push @ar,$derived->{$minor};
		}
		my $string_all_allele = join(",",@ar);
		my $popstring="c($string_all_allele)";
		return $popstring;
	}
    
	
    
}




{
    package CMHTest;
    use FindBin qw/$RealBin/;
    use lib "$RealBin/Modules";
    use Test;
    use Pileup;
    use Test::TSynchronized;
    use Test::TMaxCoverage;
    use MajorAlleles;
    

    sub runTests
    {
        run_MaxCoverageTests();
        run_SynchronizedTests();
	run_majorminor();

        exit;
    }
    
    sub run_majorminor
    {
	my($maj,$min);
	
	($maj,$min)=get_major_minor_alleles([{A=>11,T=>0,C=>0,G=>0},{A=>0,T=>10,C=>0,G=>0},{A=>0,T=>0,C=>9,G=>0},{A=>0,T=>0,C=>0,G=>9}]);
	is($maj,"A","identification of major allele; correct allele");
	is($min,"T","identification of minor allele; correct allele");
	
	($maj,$min)=get_major_minor_alleles([{A=>3,T=>0,C=>0,G=>0},{A=>3,T=>10,C=>0,G=>0},{A=>3,T=>0,C=>9,G=>0},{A=>2,T=>0,C=>0,G=>9}]);
	is($maj,"A","identification of major allele; correct allele");
	is($min,"T","identification of minor allele; correct allele");	
	
	($maj,$min)=get_major_minor_alleles([{A=>3,T=>0,C=>0,G=>0},{A=>3,T=>10,C=>0,G=>0},{A=>3,T=>2,C=>9,G=>0},{A=>2,T=>0,C=>0,G=>9}]);
	is($maj,"T","identification of major allele; correct allele");
	is($min,"A","identification of minor allele; correct allele");		
	
	($maj,$min)=get_major_minor_alleles([{A=>3,T=>0,C=>4,G=>0},{A=>3,T=>10,C=>0,G=>0},{A=>3,T=>2,C=>9,G=>0},{A=>2,T=>0,C=>0,G=>9}]);
	is($maj,"C","identification of major allele; correct allele");
	is($min,"T","identification of minor allele; correct allele");	
	    
    }

    
}



=head1 NAME

cmh-test.pl - This script calculates the Cochran-Mantel-Haenszel test for each SNP  

=head1 SYNOPSIS

 perl cmh-test.pl --input input.sync --output output.cmh --min-count 2 --min-coverage 4 --max-coverage 1000 --population 1-2,3-4,5-6 --remove-temp

=head1 OPTIONS

=over 4

=item B<--input>

The input file has to be synchronized pileup file. Mandatory parameter

=item B<--output>

The output file. Mandatory parameter

=item B<--min-count>

the minimum count of the minor allele; used for SNP identification. SNPs will be identified considering all populations simultanously. default=2

=item B<--min-coverage>

the minimum coverage; used for SNP identification, the coverage in ALL populations has to be higher or equal to this threshold, otherwise no SNP will be called. default=4

=item B<--max-coverage>

The maximum coverage; All populations are required to have coverages lower or equal than the maximum coverage; Mandatory
The maximum coverage may be provided as one of the following:

 '500' a maximum coverage of 500 will be used for all populations
 '300,400,500' a maximum coverage of 300 will be used for the first population, a maximum coverage of 400 for the second population and so on
 '2%' the 2% highest coverages will be ignored, this value is independently estimated for every population
  
=item B<--min-pvalue>

the minimum p-value cut off  to filter all snp with > min-pvalue cutoff [Optional parameter]

=item B<--population>

the pairwise comparsions which will be used for the cmh-test.
Pairwise comparisions have to be separated by a C<,> and the two populations which will be compared by a C<->. For example when the user provides 1-3,2-4 the script will compare population 1 with population 3 and population 2 with population 4;
Note also comparisions involving one population for several times are possible (e.g.: 1-2,1-3); Mandatory parameter

=item B<--remove-temp>

flag; remove the temporary files at the end; default=off

=item B<--test>

Run the unit tests for this script. 

=item B<--help>

Display help for this script

=back

=head1 Details

=head2 Input

Input is a single tab delimited file which contains a lightwight representation of every pileup file.
Every pileup file represents a population and will be parsed into a list of A-count:T-count:C-count:G-count:N-count:*-count

 2L	5002	G	0:0:0:17:0:0	0:0:0:28:0:0	0:0:0:31:0:0	0:0:0:35:0:0	0:1:0:33:0:0	0:3:0:31:0:0
 2L	5009	A	16:0:0:0:0:0	26:0:0:0:0:0	29:0:1:0:0:0	36:0:0:0:0:0	34:0:0:0:0:0	32:0:1:0:0:0
 2L	5233	G	0:0:5:46:0:0	0:0:0:43:0:0	0:0:0:60:0:0	0:0:3:61:0:0	0:0:0:56:0:0	0:0:0:48:0:0
 
 col 1: reference chromosome
 col 2: position in the reference chromosome
 col 3: reference genome base
 col 4: population 1
 col 5: population 2
 col n: population n-3
 
 population data are in the form
 A:T:C:G:N:*
 A: count of character A
 T: count of character T
 C: count of character C
 G: count of character G
 N: count of character N
 *: deletion, count of deletion
 
=head2 Output

 2L	5002	G	0:0:0:17:0:0	0:0:0:28:0:0	0:0:0:31:0:0	0:0:0:35:0:0	0:1:0:33:0:0	0:3:0:31:0:0	0.609
 2L	5009	A	16:0:0:0:0:0	26:0:0:0:0:0	29:0:1:0:0:0	36:0:0:0:0:0	34:0:0:0:0:0	32:0:1:0:0:0	0.957
 2L	5233	G	0:0:5:46:0:0	0:0:0:43:0:0	0:0:0:60:0:0	0:0:3:61:0:0	0:0:0:56:0:0	0:0:0:48:0:0	0.8088


 col 1: reference chromosome
 col 2: position in the reference chromosome
 col 3: reference genome base
 col 4: population 1
 col 5: population 2
 col n: population n-3
 col n+1: cmh p-value


=head1 Technical details

This script identifies the two major alleles for every SNPs and than runs the run Cochran mental haenszel test.
The script creates two temporary output files for R, having the extensions C<.rin> and C<.rout> which may be automatically removed using the option C<--remove-temp>.


=head2 Test statistic

You use the Cochran–Mantel–Haenszel test (which is sometimes called the Mantel–Haenszel test) for repeated tests of independence.
There are three nominal variables; you want to know whether two of the variables are independent of each other, and the third variable identifies the repeats.
For more details: http://udel.edu/~mcdonald/statcmh.html

=head1 AUTHORS

 Ram Vinay Pandey
 Robert Kofler
 Viola Nolte
 Christian Schloetterer

=cut
