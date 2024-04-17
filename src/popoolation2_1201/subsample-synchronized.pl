#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Path;
use File::Basename; # to get the file path, file name and file extension
use FindBin qw/$RealBin/;
use lib "$RealBin/Modules";
use List::Util qw[min max];
use MaxCoverage;
use Subsample;
use Synchronized;
use SynchronizeUtility;

my $input;
my $output="";
my $help=0;
my $test=0;
my $targetcoverage=0;
my $usermaxcoverage=0; # the user may provide one of the following: 500 or 500,400,300 or 2%
my $method; #withreplace, withoutreplace, fraction


GetOptions(
    "input=s"	        =>\$input,
    "output=s"          =>\$output,
    "max-coverage=s"    =>\$usermaxcoverage,
    "target-coverage=i" =>\$targetcoverage,
    "method=s"          =>\$method,
    "test"              =>\$test,
    "help"	        =>\$help
) or pod2usage(-msg=>"Wrong options",-verbose=>1);

# too many arguments should result in an error
pod2usage(-msg=>"Wrong options",-verbose=>1) if @ARGV;
pod2usage(-verbose=>2) if $help;
SubsampleTests::runTests() if $test;
pod2usage(-msg=>"Please provide an existing input file",-verbose=>1) unless -e $input;
pod2usage(-msg=>"Please provide an output file",-verbose=>1) unless $output;
pod2usage(-msg=>"Please provide a maximum coverage",-verbose=>1) unless $usermaxcoverage;
pod2usage(-msg=>"Please provide a target coverage",-verbose=>1) unless $targetcoverage;
pod2usage(-msg=>"Please provide a sampling method",-verbose=>1) unless $method;


my $paramfile=$output.".params";
open my $pfh, ">",$paramfile or die "Could not open $paramfile\n";
print $pfh "Using input\t$input\n";
print $pfh "Using output\t$output\n";
print $pfh "Using maximum coverage\t$usermaxcoverage\n";
print $pfh "Using subsample method\t$method\n";
print $pfh "Using target coverage\t$targetcoverage\n";
print $pfh "Using test\t$test\n";
print $pfh "Using help\t$help\n";
close $pfh;

my $maxcoverage=get_max_coverage($input,$usermaxcoverage);
my $pp=get_sumsnp_synparser(1,$targetcoverage,$maxcoverage);
my $subsampler=get_subsampler($method,$targetcoverage);

open my $ifh, "<",$input or die "Could not open input file $input";
open my $ofh, ">",$output or die "Could not open output file $output";

while(my $line=<$ifh>)
{
    chomp $line;
    next unless $line;
    my $p=$pp->($line);
    next unless $p->{iscov};
    
    my $samplecount=@{$p->{samples}};
    for(my $i=0; $i<$samplecount; $i++)
    {
        # randomly subsample every sample to the given coverage        
        my $subsampled=$subsampler->($p->{samples}[$i]);
        $p->{samples}[$i]=$subsampled;
    }
    my $toprint=format_synchronized($p);
    print $ofh $toprint."\n";
}

{
    use strict;
    use warnings;
    package SubsampleTests;
    use FindBin qw/$RealBin/;
    use lib "$RealBin/Modules";
    use Test::TSubsample;
    use Test::TSynchronized;
    use Test::TMaxCoverage;
    
    
    sub runTests
    {
        run_MaxCoverageTests();
        run_SynchronizedTests();
        run_SubsampleTests();
        exit;
    }
}



=head1 NAME

perl subsample-synchronized.pl - Reduce the coverage of a synchronized file, by random subsampling, to the given target coverage

=head1 SYNOPSIS

perl subsample-synchronized.pl --input input.sync --output output.sync --target-coverage 50 --max-coverage 2%  --method withoutreplacement

=head1 OPTIONS

=over 4

=item B<--input>

The input file in the synchronized format; Mandatory.

=item B<--output>

The output file, will be a synchronized file  Mandatory.

=item B<--target-coverage>

Reduce the coverage of the pileup-file to the here provided value; The target coverage also acts as minimum coverage,
i.e.: if the coverage in any population is smaller than the targetcoverage the whole pileup entry is discarded. Mandatory

=item B<--max-coverage>

The maximum coverage; All populations are required to have coverages lower or equal than the maximum coverage; Mandatory
The maximum coverage may be provided as one of the following:

 '500' a maximum coverage of 500 will be used for all populations
 '300,400,500' a maximum coverage of 300 will be used for the first population, a maximum coverage of 400 for the second population and so on
 '2%' the 2% highest coverages will be ignored, this value is independently estimated for every population

=item B<--method>

Specify the method for subsampling of the synchronized file. Either: withreplace, withoutreplace, fraction; Mandatory
 
 withreplace: subsample with replacement
 withoutreplace: subsample without replacement
 fraction: calculate the exact fraction of the allele frequencies and linearly scale them to the C<--target-coverage> and rounding to an integer; 

=item B<--help>

Display help for this script

=back

=head1 Details

=head2 Input synchronized

A synchronized file

=head2 Output

The output will be a reduced coverage synchronized file

=cut



