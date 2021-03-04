#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Path;
use File::Basename; # to get the file path, file name and file extension
use FindBin qw/$RealBin/;
use lib "$RealBin/Modules";



my $input="";
my $output="";
my $windowsize=1;
my $logtransform=0;
my $test=0;
my $help=0;

GetOptions(
    "input=s"	        =>\$input,
    "output=s"          =>\$output,
    "window-size=i"     =>\$windowsize,
    "log-transform"     =>\$logtransform,
    "test"              =>\$test,
    "help"	        =>\$help
    );

pod2usage(-verbose=>2) if $help;
die "Tests not implemented for this script" if $test;
pod2usage(-msg=>"Input file does not exist",-verbose=>1) unless -e $input;
pod2usage(-msg=>"No output file has been provided",-verbose=>1) unless $output;
pod2usage(-msg=>"Output file has to have .igv extension!",-verbose=>1) unless $output=~m/\.igv$/;



my($header,$reader)=get_pwc_reader($input);
open my $ofh,">",$output or die "Could not open output file";

print $ofh format_header($header)."\n";

while(my $e=$reader->())
{
    
    $e=set_start_end($e,$windowsize);
    $e=log_transform($e) if $logtransform;
    my $toprint=format_data($e);
    print $ofh $toprint."\n";
}


exit;

{
    use strict;
    use warnings;
    
    sub set_start_end
    {
        my $e=shift;
        my $windowsize=shift;
        # igv is zero-based
        # pileup is one-based
        my $pos=$e->{pos}-1;
        my $start=$pos-($windowsize/2);
        my $end=$pos+($windowsize/2);
        $start  = _intround($start);
        $end    = _intround($end);
        $e->{start}=$start;
        $e->{end}=$end;
        return $e;
    }
    
    sub format_data
    {
        my $e=shift;
        my @a=();
        push @a,$e->{chr};
        push @a,$e->{start};
        push @a,$e->{end};
        push @a,"snp";
        foreach my $t (@{$e->{samples}})
        {
            push @a,$t;
        }
        return join("\t",@a);
    }
    
    sub format_header
    {
        my $header=shift;
        #Chromosome 	Start 	End 	Feature 	Patient-One 	Patient-Two 	Patient-Three
        my @toprint=();
        push @toprint,"Chromosome";
        push @toprint,"Start";
        push @toprint,"End";
        push @toprint,"Feature";
        
        foreach my $t (@$header)
        {
            push @toprint,$t;
        }
        
        return join("\t",@toprint);
    }
    
    sub log_transform
    {
        my $e=shift;
        my $samples=$e->{samples};
        
        for (my $i=0;$i<@$samples; $i++)
        {
            next if $samples->[$i]==0;
            $samples->[$i]=-1*log($samples->[$i])/log(10);
        }
        return $e;
    }
    
    sub get_pwc_reader
    {
        my $file=shift;
        
        # get the header information
    # 2R	2358	1	1.000	88.0	1:2=0.00434153	1:3=0.00000000	1:4=0.00434153	2:3=0.00434153	2:4=0.00000000	3:4=0.00434153
    # 2R	2364	1	1.000	94.0	1:2=0.00180180	1:3=0.00000000	1:4=0.00180180	2:3=0.00180180	2:4=0.00000000	3:4=0.00180180
    # 2R	2367	1	1.000	95.0	1:2=0.00178253	1:3=0.00000000	1:4=0.00178253	2:3=0.00178253	2:4=0.00000000	3:4=0.00178253
    # 2R	2369	1	1.000	95.0	1:2=0.00401606	1:3=0.00000000	1:4=0.00401606	2:3=0.00401606	2:4=0.00000000	3:4=0.00401606
        open my $ifh,"<",$file or die "Could not open input file";
        my $line=<$ifh>;
        close $ifh;
        chomp $line;
        my @a=split /\t/,$line;
        shift @a; shift @a; shift @a; shift @a; shift @a;
        
        my $header=[];
        foreach my $e(@a)
        {
            my($th)=$e=~m/^([^=]+)/;
            push @$header,$th;
        }
        
        open $ifh,"<",$file or die "Could not open input file";
        return ($header,sub
        {
            my $line=<$ifh>;
            return undef unless $line;
            chomp $line;
            my @a=split /\t/,$line;
            my $chr=shift @a;
            my $pos=shift @a;
            my $snps=shift @a;
            my $coveredfraction=shift @a;
            my $mincovered=shift @a;
            my $samples=[];
            foreach my $e (@a)
            {
                my($temp)=$e=~m/=(.+)/;
                push @$samples,$temp;
            }
            
            return {
                chr=>$chr,
                pos=>$pos,
                snpcount=>$snps,
                coveredfraction=>$coveredfraction,
                mincovered=>$mincovered,
                samples=>$samples
            };
        })
    }
    
    sub _intround
    {
        
        my $i=shift;            #2.5
        my $integer=int($i);    #2
        my $rest=$i-$integer;   #0.5
        if($rest<0.5)
        {
            return $integer;    #2
        }
        else{
            return $integer+1;  #3
        }
    }
}

=head1 NAME

pwc2igv - convert the results of pairwise comparisions (Fst, fisher-test) into the igv-format (for display in IGV)

=head1 SYNOPSIS

 pwc2igv.pl --input fisher.results --output fisher.results.igv --log-transform 

=head1 OPTIONS

=over 4

=item B<--input>

The input file. Has to be the output of either C<fst-sliding.pl> or C<fisher-test.pl>. Mandatory parameter

=item B<--output>

The output file. Must have the extension C<.igv>. Can be directly loaded into the IGV.  Mandatory parameter

=item B<--window-size>

specify the size of the window; optional parameter; default=1

=item B<--log-transform>

flag; log transform the measure (-log10). This may be especially useful for the results of <fisher-test.pl>; optional parameter; default=off

=item B<--test>

Run the unit tests for this script. 

=item B<--help>

Display help for this script

=back

=head1 DETAILS

=head2 INPUT

The input needs to be a file created either by C<fisher-test.pl> or by C<--fst-sliding>, e.g.:

 2L      68500   360     1.000   62.1    1:2=0.01873725  1:3=0.02131245  2:3=0.01521177
 2L      69500   118     1.000   71.9    1:2=0.00969479  1:3=0.00116059  2:3=0.00905794
 2L      70500   269     1.000   63.6    1:2=0.01955417  1:3=0.01547995  2:3=0.01300569

 col1: reference contig (chromosome)
 col2: mean position of the sliding window
 col3: number of SNPs found in the window (not considering sites with a deletion) 
 col4: fraction of the window which has a sufficient coverage (min. coverage <= cov <= max. coverage) in every population;
 col5: average minimum coverage in all populations

 ....
 1:2 the pairwise Fst/fisher-pvalue for population 1 and 2
 1:3 the pairwise Fst/fisher-pvalue for population 1 and 3

=head2 OUTPUT

An .igv file, for details see: http://www.broadinstitute.org/software/igv/IGV

=head1 AUTHORS

Robert Kofler

Christian Schloetterer

=cut