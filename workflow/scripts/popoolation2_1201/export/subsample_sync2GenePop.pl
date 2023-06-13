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
my $mincount=1;
my $diploid=0;
my $userregion="";
my $usermaxcoverage=0; # the user may provide one of the following: 500 or 500,400,300 or 2%
my $method; #withreplace, withoutreplace, fraction

# --input /Users/robertkofler/dev/testfiles/sync100000.sync --output /Users/robertkofler/dev/testfiles/output/pseudofasta --method withoutreplace --target-coverage 10 --max-coverage 400 --region 2L:10000-10200

GetOptions(
    "input=s"	        =>\$input,
    "output=s"          =>\$output,
    "max-coverage=s"    =>\$usermaxcoverage,
    "target-coverage=i" =>\$targetcoverage,
    "min-count=i"       =>\$mincount,
    "method=s"          =>\$method,
    "region=s"          =>\$userregion,
    "diploid"           =>\$diploid,
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


my $paramfile=$output.".params";
open my $pfh, ">",$paramfile or die "Could not open $paramfile\n";
print $pfh "Using input\t$input\n";
print $pfh "Using output\t$output\n";
print $pfh "Using minimum count\t$mincount\n";
print $pfh "Using maximum coverage\t$usermaxcoverage\n";
print $pfh "Using subsample method\t$method\n";
print $pfh "Using target coverage\t$targetcoverage\n";
print $pfh "Using region\t$userregion\n";
print $pfh "Using diploid data\t$diploid\n";
print $pfh "Using test\t$test\n";
print $pfh "Using help\t$help\n";
close $pfh;

my $maxcoverage=get_max_coverage($input,$usermaxcoverage);
my $pp=get_sumsnp_synparser($mincount,$targetcoverage,$maxcoverage);
my $subsampler=get_subsampler($method,$targetcoverage);
my $popcount=get_popcount_forsyncfile($input);

my($regionchr,$regionstart,$regionend,$region)=Utility::parse_region($userregion);

open my $ifh, "<",$input or die "Could not open input file $input";
my $activeflag=0;
my $snpsinregion=[];
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
    next unless $p->{ispuresnp};
    my $samplecount=@{$p->{samples}};
    die "Numbers of populations is not consistent within the file  $samplecount vs $popcount" unless $samplecount eq $popcount;
    push @$snpsinregion,$p;
}

# reorganize data and convert them into GenPop genotypes
my $genepopar=[];
for my $i (0..($popcount-1))
{
    my $popar=[];
    foreach my $snp(@$snpsinregion)
    {
        my $pop=$snp->{samples}[$i];
        my $subs=$subsampler->($pop);
        my $gp=Utility::convert_to_genepopencoding($subs,$diploid);
        push @$popar,$gp;
    }
    die "User specified region does not contain any SNPs, can not create GenePop output for no SNPs" unless @$popar;
    push @$genepopar,$popar;
}


Utility::print_genepop($output,$genepopar,$userregion,$targetcoverage);

    #for(my $i=0; $i<$popcount; $i++)
    #{
    #    # randomly subsample every sample to the given coverage        
    #    my $subsampled=$subsampler->($p->{samples}[$i]);
    #
    #}
close $ifh;




exit;

{
    use strict;
    use warnings;
    package Utility;
    
    sub print_genepop
    {
        my $output=shift;
        my $gpa=shift;
        my $region=shift;
        my $coverage=shift;
        
        my $endl="\n";
        open my $ofh, ">", $output or die "Could not open output file";
        
        my $snpcount=@{$gpa->[0]};
        my $toprint=[];
        foreach my $i (1..$snpcount)
        {
            push @$toprint,"Loc$i";
        }
        print $ofh "Haploid genotypes for region $region, subsampled to a coverage of $coverage; NOTE: haplotype structure is a mere artifact!".$endl;
        print $ofh join(", ",@$toprint).$endl;
        
        foreach my $pop (@$gpa)
        {
            print $ofh "Pop".$endl;
            INDIVIDUALS: while(1)
            {
                $toprint=["gt ,"];
                foreach my $snp (@$pop)
                {
                    last INDIVIDUALS unless @$snp;
                    my $gt=pop @$snp;
                    push @$toprint,$gt;
                }
                print $ofh join(" ",@$toprint).$endl;
            }
        }
    }
    
    sub convert_to_genepopencoding{
        my $e=shift;
        my $diploid=shift;
        
        my $cov=$e->{A}+$e->{T}+$e->{C}+$e->{G}+$e->{N}+$e->{del};
        my $toret;
        if($diploid)
        {
            push @$toret,"0101" foreach (1..$e->{A});
            push @$toret,"0202" foreach (1..$e->{T});
            push @$toret,"0303" foreach (1..$e->{C});
            push @$toret,"0404" foreach (1..$e->{G});
            push @$toret,"0000" foreach (1..$e->{N});
            push @$toret,"0000" foreach (1..$e->{del});
        }
        else
        {
            push @$toret,"01" foreach (1..$e->{A});
            push @$toret,"02" foreach (1..$e->{T});
            push @$toret,"03" foreach (1..$e->{C});
            push @$toret,"04" foreach (1..$e->{G});
            push @$toret,"00" foreach (1..$e->{N});
            push @$toret,"00" foreach (1..$e->{del});
        }

        
        die "improper size"unless @$toret == $cov;
        return $toret;
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

perl subsample-sync2GenePop.pl - Convert a given region from a synchronized file into a GenePop file (may be converted into files compatible with Arlequin)

=head1 SYNOPSIS

perl subsample-sync2GenePop.pl --input input.sync --output output.genepop --target-coverage 20 --max-coverage 2%  --method withoutreplace --region 2R:10000-20000,40000-50000 --min-count 2

=head1 OPTIONS

=over 4

=item B<--input>

The input file in the synchronized format; Mandatory.

=item B<--output>

The output file, will be a GenePop file;  Mandatory.
NOTE: haplotype information is lost during pooling, therefore the haplotypes found in the GenePop file are artifactual!

=item B<--min-count>

This script only works with SNPs. To identify a SNP a minimum allele count is required. default=1

=item B<--target-coverage>

Reduce the coverage per population to this value; The target coverage also acts as minimum coverage,
i.e.: if the coverage in any population is smaller than the targetcoverage the whole entry is discarded.
The targetcoverage will also be the number of samples per population in the GenePop file!
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

=item B<--diploid>

Flag; switch to diploid encoding of the GenePop file; per default it is producing haploid encoding. This option may be important for some convertes as they do not recognize the haploid GenePop file. default=off

=item B<--help>

Display help for this script

=back

=head1 Details

=head2 Input synchronized

A synchronized file

=head2 Output

A file compatible with GenePop http://genepop.curtin.edu.au/index.html
The GenePop webserver also allows conversion into files compatible with: Arelequin, Fstat, Linkdos, Biosys

The script produces haploid genotypes as pooling does not allow infering diploid genotypes:
Following a slightly modified example from the Genepop web page:

 Microsatelite loci on Chiracus radioactivus, a pest species
 Loc1, Loc2, Loc3
 Pop
 AA8, 04 07 03
 AA9, 04 06 02
 A10, 02 06 01
 A11, 04 06 01
 Pop
 AF, 00 00 00
 AF, 02 03 01
 AF, 02 03 02
 AF, 02 03 00

=head2 Disclaimer

NOTE: during pooling and sequencing the haplotype information is lost. Thus any haplotype structure found in the GenePop file is merley an artefact!
Do not use this GenePop file for any analysis involving haplotypes. However, it may be used to measure differentiation between populations/subpopualtions or to calculate Tajima's D, Pi etc

=cut



