#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

my $input;
my $output;
my $help=0;

my $minpvalue=1.0e-20;


GetOptions(
    "input=s"	    =>\$input,
    "output=s"	    =>\$output,
    "min-pvalue=s"  =>\$minpvalue,
    "help"	    =>\$help
) or pod2usage(-msg=>"Wrong options",-verbose=>1);

pod2usage(-verbose=>2) if $help;
pod2usage(-msg=>"A input file has to be provided\n",-verbose=>1) unless -e $input;
pod2usage(-msg=>"A output file has to be provided\n",-verbose=>1) unless $output;
pod2usage(-msg=>"Output needs file extension .gwas\n",-verbose=>1) unless $output=~/\.gwas$/;

open my $ifh, "<", $input or die "Could not open input file";
open my $ofh, ">", $output or die "Could not open output file";

my $counter=1;
print $ofh "CHR\tBP\tSNP\tP\n";
while(my $line=<$ifh>)
{
    chomp $line;
    
    my @ar=split /\t/,$line;
    my $chr=shift @ar;
    my $pos=shift @ar;
    my $pvalue=pop @ar;
    $pvalue=$minpvalue if $pvalue< $minpvalue;
    print $ofh "$chr\t$pos\tsnp$counter\t$pvalue\n";
    $counter++;
}
close $ofh;
close $ifh;
    #CHR: chromosome (aliases chr, chromosome)
    #BP: nucleotide location (aliases bp, pos, position)
    #SNP: SNP identifier (aliases snp, rs, rsid, rsnum, id, marker, markername)
    #P: p-value for the association (aliases p, pval, p-value, pvalue, p.value)


#2R	4459	N	0:55:45:0:0:0	0:40:54:2:0:0	0:55:45:0:0:0	0:40:54:2:0:0	0.01910801
#2R	9728	N	0:56:44:0:0:0	0:44:55:0:0:0	0:56:44:0:0:0	0:44:55:0:0:0	0.0278406
#2R	9828	N	0:56:44:0:0:0	1:43:56:0:0:0	0:56:44:0:0:0	1:43:56:0:0:0	0.01637122
#2R	9928	N	44:0:1:55:0:0	56:0:0:43:0:0	44:0:1:55:0:0	56:0:0:43:0:0	0.02111846

=head1 NAME

cmh2gwas.pl - Converts the output of cmh-test.pl into a gwas file format accepted by the IGV

=head1 SYNOPSIS

 perl cmh2gwas.pl --input input.cmh --output output.gwas

=head1 OPTIONS

=over 4

=item B<--input>

The input file has to be synchronized pileup file. Mandatory parameter

=item B<--output>

The output file. Mandatory parameter

=item B<--min-pvalue>

IGV has some problems displaying very low p-values; Using this option all p-values being smaller than C<--min-pvalue> will be set to C<--min-pvalue>; default=1.0e-20 

=item B<--help>

Display help for this script

=back

=head1 Details

=head2 Input

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
 
=head2 Output






=head1 AUTHORS

 Robert Kofler
 Christian Schloetterer

=cut
