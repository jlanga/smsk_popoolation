#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($RealBin);
use lib "$RealBin/Modules";
use SynchronizeUtility;
use Pileup;


my $input="";
my $output="";
my $fastqtype="illumina";
my $minQual=20;
my $test=0;
my $help=0;

GetOptions(
    "input=s"           =>\$input,
    "output=s"          =>\$output,
    "fastq-type=s"      =>\$fastqtype,
    "min-qual=i"        =>\$minQual,
    "test"              =>\$test,
    "help"              =>\$help
) or die "Invalid arguments";

pod2usage(-verbose=>2) if $help;
VarTest::runTests() if $test;
pod2usage(-msg=>"Could not find input file",-verbose=>1) unless -e $input;
pod2usage(-msg=>"Output file not provided",-verbose=>1) unless  $output;
pod2usage(-msg=>"Mininum quality needs to be >0",-verbose=>1) unless $minQual;


# Writing the .param file
my $paramfile=$output.".params";
open my $pfh, ">",$paramfile or die "Could not open $paramfile\n";
print $pfh "Using input\t$input\n";
print $pfh "Using output\t$output\n";
print $pfh "Using fastq-type\t$fastqtype\n";
print $pfh "Using min-qual\t$minQual\n";
print $pfh "Using help\t$help\n";
print $pfh "Using test\t$test\n";
close $pfh;

my $pp=get_basic_mpileupparser($fastqtype,$minQual);

open my $ifh, "<", $input or die "Could not open input file $input";
open my $ofh, ">", $output or die "Could not open output file $output";


while(my $l=<$ifh>)
{
    chomp $l;
    my $p=$pp->($l);
    my @ar=();
    my $entries=$p->{entries};
    push @ar,$p->{chr};
    push @ar,$p->{pos};
    push @ar,$p->{refc};
    foreach my $e (@$entries)
    {
        my $f=format_parsed_pileup($e);
        push @ar,$f;
    }
    my $toprint =join("\t",@ar);
    print $ofh $toprint."\n";
}


{
    package VarTest;
    use FindBin qw($RealBin);
    use lib "$RealBin/Modules";
    use Test::PileupParser;
    use Pileup;
    use Test;
    sub runTests
    {
        run_PileupParserTests();
        exit;
    }
    
}


=head1 NAME

mpileup2sync.pl - Converts a mpileup (multi-pileup) file into the synchronized file format

=head1 SYNOPSIS

 mpileup2sync.pl --input severalpops.mpileup --min-qual 20 --output output.txt

=head1 OPTIONS

=over 4

=item B<--input>

The input file in mpileup format. At least two have to be specified. Mandatory parameter

=item B<--output>

The output file. Mandatory parameter

=item B<--min-qual>

the minimum quality of a base to be reported; default=20

=item B<--fastq-type>

The encoding of the quality characters; Must either be 'sanger' or 'illumina'; 

 Using the notation suggested by Cock et al (2009) the following applies:
 'sanger'   = fastq-sanger: phred encoding; offset of 33
 'solexa'   = fastq-solexa: -> NOT SUPPORTED
 'illumina' = fastq-illumina: phred encoding: offset of 64
 
 See also:
 Cock et al (2009) The Sanger FASTQ file format for sequecnes with quality socres,
 and the Solexa/Illumina FASTQ variants; 

default=illumina

=item B<--test>

Run the unit tests for this script. 

=item B<--help>

Display help for this script

=back

=head1 DETAILS

Combine different pileup files

=head2 INPUT

a mpileup file as produced by samtools, for example:
 
 YHet    4067    N       9       ttttTtttt       aaab_Za_b       2       tt      `b      2       tt      \b
 YHet    4068    N       9       c$cccCcccc      a_a_a]a_`       2       cc      ab      2       cc      `b
 YHet    4069    N       8       aaaAaaaa        ]aaa_^__        2       aa      ab      2       aa      \b
 YHet    4070    N       8       a$a$aAaaaa      \a]^YX_a        2       a$a     ab      2       aa      Wb

=head2 OUTPUT

Output is a single tab delimited file which contains a lightwight representation of every mpileup file.
Every pileup file represents a population and will be parsed into a list of A-count:T-count:C-count:G-count:N-count:*-count
The order of samples as in the mpileup file will be preserved, i.e.: the first population in the mpileup file will be the first population in the synchronized file and so on. 

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

=head1 AUTHORS

Robert Kofler

Christian Schloetterer

=cut