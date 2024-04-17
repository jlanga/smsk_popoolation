#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Path;
use File::Basename; # to get the file path, file name and file extension
use FindBin qw/$RealBin/;
use lib "$RealBin/../Modules";
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
my $userregion="";
my $usermaxcoverage=0; # the user may provide one of the following: 500 or 500,400,300 or 2%
my $method; #withreplace, withoutreplace, fraction
my $regionid="";

# --input /Users/robertkofler/dev/testfiles/sync100000.sync --output /Users/robertkofler/dev/testfiles/output/pseudofasta --method withoutreplace --target-coverage 10 --max-coverage 400 --region 2L:10000-10200

GetOptions(
    "input=s"	        =>\$input,
    "output=s"          =>\$output,
    "max-coverage=s"    =>\$usermaxcoverage,
    "target-coverage=i" =>\$targetcoverage,
    "method=s"          =>\$method,
    "region=s"          =>\$userregion,
    "region-id=s"       =>\$regionid,
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
pod2usage(-msg=>"Please provide a region which should be converted to a multiple fasta",-verbose=>1) unless $userregion; 
$regionid=$userregion unless $regionid;


my $paramfile=$output.".params";
open my $pfh, ">",$paramfile or die "Could not open $paramfile\n";
print $pfh "Using input\t$input\n";
print $pfh "Using output\t$output\n";
print $pfh "Using maximum coverage\t$usermaxcoverage\n";
print $pfh "Using subsample method\t$method\n";
print $pfh "Using target coverage\t$targetcoverage\n";
print $pfh "Using region\t$userregion\n";
print $pfh "Using region-id\t$regionid\n";
print $pfh "Using test\t$test\n";
print $pfh "Using help\t$help\n";
close $pfh;

my $maxcoverage=get_max_coverage($input,$usermaxcoverage);
my $pp=get_sumsnp_synparser(1,$targetcoverage,$maxcoverage);
my $subsampler=get_subsampler($method,$targetcoverage);
my $popcount=get_popcount_forsyncfile($input);

my($regionchr,$regionstart,$regionend,$region)=Utility::parse_region($userregion);

open my $ifh, "<",$input or die "Could not open input file $input";
my $activeflag=0;
SYNCLINE: while(my $line=<$ifh>)
{
    chomp $line;
    next unless $line;
    my($chr,$pos)=split /\s+/,$line;
    
    # set activeflag
    if($activeflag)
    {
        
        last SYNCLINE if($pos >$regionend or $chr ne $regionchr)
    }
    unless($activeflag)
    {
        $activeflag=1 if($chr eq $regionchr  and $pos >=$regionstart);
    }
    #skip everthing which is not in the user-defined region
    next unless $activeflag;
    next unless(exists($region->{$chr}{$pos}));
    die "Unallowed state chromosome not equal to targetchromosome $chr vs $regionchr" unless $chr eq $regionchr;
    

    my $p=$pp->($line);
    next unless $p->{iscov};
    my $samplecount=@{$p->{samples}};
    die "Numbers of populations is not consistent within the file  $samplecount vs $popcount" unless $samplecount eq $popcount;
    
    my $tostr="";
    for(my $i=0; $i<$popcount; $i++)
    {
        # randomly subsample every sample to the given coverage        
        my $subsampled=$subsampler->($p->{samples}[$i]);
        $tostr.=syncsample2string($subsampled);
    }
    
    my @strar=split //,$tostr;
    
    $region->{$chr}{$pos}=\@strar;
}
close $ifh;

## fill in the missing values;
# missing due to coverage issues
my @nfiller= split //, "N" x ($popcount * $targetcoverage);
for my $i ($regionstart..$regionend)
{
    next unless(exists($region->{$regionchr}{$i}));
    my $count= @{$region->{$regionchr}{$i}};
    next if $count;
    my @temp=@nfiller;
    $region->{$regionchr}{$i}=\@temp;
}


## print one sequence after the other

my $counter=0;
open my $ofh, ">", $output or die "Could not open output file";
FASTA: while(1)
{
    my $fasta=[];
    my $popcount=int($counter/$targetcoverage);
    my $samplecount=$counter % $targetcoverage;
    $samplecount++;
    $popcount++;
    
    for my $i ($regionstart..$regionend)
    {
        
        next unless(exists($region->{$regionchr}{$i}));
        my $tar=$region->{$regionchr}{$i};
        unless(scalar(@$tar))
        {
            last FASTA;
        }
        my $char=shift @$tar;
        push @$fasta,$char;
    }
    
    Utility::write_fasta_entry($ofh,$regionid,$popcount,$samplecount,$fasta);
    $counter++;  
}





exit;

{
    use strict;
    use warnings;
    package Utility;
    
    sub write_fasta_entry{
        my $ofh=shift;
        my $regionid=shift;
        my $popcount=shift;
        my $samplecount=shift;
        my $fasta=shift;
        my $topr=join("",@$fasta);
        
        print $ofh ">$regionid"."_pop$popcount"."_sample$samplecount\n";
        print $ofh $topr."\n";
    }
    
    sub parse_region
    {
        my $region=shift; #2R:123-145
        
        die "Region invalid $region; must be of the form chr:start-end,start-end" unless $region=~m/:/;
        die "Region invalid $region; must be of the form chr:start-end,start-end" unless $region=~m/[-]/;
        my $lowest=undef;
        my $highest=undef;
        my $poscollection={};
        
        my($chr,$temp)=split /:/,$region;
        my @ar=split/,/,$temp;
        push @ar,$temp unless(@ar);# if only a single exon was provided
        
        
        
        foreach my $a (@ar)
        {
            die "Region invalid $region; must be of the form chr:start-end,start-end"unless $region=~m/[-]/;
            my($start,$end)=split /-/,$a;
            die "Region invalid: $region Start position must be smaller than the end position" if $start > $end;
            $lowest =$start if(not defined($lowest));
            $highest=$end if(not defined($highest));
            $lowest =$start if $start <$lowest;
            $highest=$end if $end > $highest;
            
            for my $i($start..$end)
            {
                $poscollection->{$chr}{$i}=[];
            }
            
        }
        return ($chr,$lowest,$highest,$poscollection);
        
    }
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

perl subsample-sync2fasta.pl - Convert a given region from a synchronized file into a multi-fasta file for use with third party software

=head1 SYNOPSIS

perl subsample-sync2fasta.pl --input input.sync --output output.sync --target-coverage 20 --max-coverage 2%  --method withoutreplacement --region 2R:10000-20000,40000-50000 --region-id supercoolgene

=head1 OPTIONS

=over 4

=item B<--input>

The input file in the synchronized format; Mandatory.

=item B<--output>

The output file, will be a multi-fasta file;  Mandatory.
NOTE: haplotype information is lost during pooling, therefore the haplotypes found in the fasta file are artifactual!

=item B<--target-coverage>

Reduce the coverage per population to this value; The target coverage also acts as minimum coverage,
i.e.: if the coverage in any population is smaller than the targetcoverage the whole entry is discarded.
The targetcoverage will also be the number of samples per population in the fasta file!
Mandatory

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

=item B<--region>

Specify the region which should be subsampled; Has to be in the format: 'chr:start-end,start-end,start-end' Mandatory
For example: '2R:1-100,200-220,240-280' NOTE that the option to provide several regions separated by a comma ',' allows to specify exonic regions

=item B<--region-id>

Provide a short and descriptive id of the requested region, for example a gene name; default=C<--region>

=item B<--help>

Display help for this script

=back

=head1 Details

=head2 Input synchronized

A synchronized file

=head2 Output

The output will be a multiple fasta file which could be used with third party software for Population Genetics.
The output file will contain information about the populations and the number of samples used.
In the following example the synchronized file contained two populations and the C<--target-coverage> (thus also the number of samples) was set to 2:

 >2L:10000-10020_pop1_sample1
 GGAGGAGAATGCAAAAAAGCT
 >2L:10000-10020_pop1_sample2
 GGAGGAGAATGCAAAAAAGCT
 >2L:10000-10020_pop2_sample1
 GGAGGAGAATGCAAAAAAGCT
 >2L:10000-10020_pop2_sample2
 GGAGGAGAATGCAAAAAAGCT

=head2 Disclaimer

NOTE: during pooling and sequencing the haplotype information is lost. Thus any haplotype structure found in the fasta file is merley an artifact!
Do not use this multifasta file for any analysis involving haplotypes. However, it may be used to measure differentiation between populations/subpopualtions or to calculate Tajima's D, Pi etc

=cut



