#!/usr/bin/env perl
#
# Authors:
# Dr. Robert Kofler
# Ram Vinay Panday
# 02.02.2010
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Path;
use File::Basename; # to get the file path, file name and file extension
use List::Util qw[min max];
use POSIX qw(ceil floor);
use FindBin qw($RealBin);
use lib "$RealBin/Modules";
use SynchronizeUtility;
use Pileup;


#--input test/t.pile --input test/r.pile --input test/s.pile --output test/test.comb --fastq-type illumina

my @input = ();
my $minQual=20;
my $outputFile="";
my $help=0;
my $test=0;
my $reportat=100000;
my $reportcounter=0;
my $fastqtype="illumina";

GetOptions(
    "input=s"	    =>\@input,
    "min-qual=s"    =>\$minQual,
    "output=s"	    =>\$outputFile,
    "fastq-type=s"  =>\$fastqtype,
    "test"          =>\$test,
    "help"	    =>\$help
) or pod2usage(-msg=>"Wrong options",-verbose=>1);

pod2usage(-verbose=>2) if $help;
SyncTest::runTests() if $test;
pod2usage(-msg=>"Not enough input files, have to be two or more",-verbose=>1) if @input<2;
foreach my $ifile (@input)
{
    pod2usage(-msg=>"Input file does not exist",-verbose=>1) unless -e $ifile;    
}
pod2usage(-msg=>"No output file has been provided",-verbose=>1) unless $outputFile;
$fastqtype=lc($fastqtype);


my $paramfile=$outputFile.".params";
open my $pfh, ">",$paramfile or die "Could not open $paramfile\n";
print $pfh "Using input\t$_\n" foreach @input;
print $pfh "Using min-qual\t$minQual\n";
print $pfh "Using output\t$outputFile\n";
print $pfh "Using fastq-type\t$fastqtype\n";
print $pfh "Using test\t$test\n";
print $pfh "Using help\t$help\n";
close $pfh;


# qualencoding,mincount,mincov,maxcov,minqual
my $pp=get_pileup_parser($fastqtype,1,1,2000000000,$minQual);


open my $ofh, ">$outputFile" or die "Could not open output file";

my $filecount=@input;
my $pileupFiles;

# create an instance of a Pileup file for each
push @$pileupFiles,PileupFile->new($_,$pp) foreach @input;

# create a synchronizer from the pileup files
my $synchronizer=Synchronizer->new($pileupFiles);


while(my $entry=$synchronizer->next())


{
    # print the lines which  are far away
    Utility::printLine($ofh,$entry,$filecount);
    
    # increas the counter
    $reportcounter++;
    print "Processed $reportcounter positions\n" unless ($reportcounter % $reportat);
    
}


exit;



{
    package Utility;
    use strict;
    use warnings;
    use FindBin qw($RealBin);
    use lib "$RealBin/Modules";
    use SynchronizeUtility;
    
    
    sub printLine
    {
        my $ofh=shift;
        my $entry=shift;
        my $filecount=shift;
        
        my $h;
        foreach(@{$entry->{data}})
        {
            $h->{$_->{index}}=$_;
        }
        
        
        print $ofh "$entry->{chr}\t$entry->{pos}\t$entry->{refc}";
        
        for(my $i=0; $i<$filecount; $i++)
        {
            
            my $e=$h->{$i};
            my $str="\t".format_parsed_pileup($e);
            print $ofh $str;
        }
        print $ofh "\n";
    }
}

{
    package Synchronizer;
    use strict;
    use warnings;
    
    sub new
    {
        my $self=shift;
        my $pileupfiles=shift;
        
        my $curLines=[];
        
        # initialize every pileup file
        foreach (@$pileupfiles)
        {
            my $e=$_->nextLine();
            push @$curLines,$e if $e; # skip if empty
        }
        
        my $e= {
            pileupfiles=>$pileupfiles,
            curlines=>$curLines
        };
        
        return bless $e,  __PACKAGE__;
        
    }
    
    
    sub next
    {
        my $self=shift;
        
        my($overhead,$toprint)=$self->_splitCurrentLines();
        
        
        # read from the files
        my $novel=$self->_getNewLines($toprint);
        
        # add the newly read files to the reminder
        push @$overhead,@$novel;
        $self->{curlines}=$overhead;
        
        die "toprint emtpy while overhead still contains entries" if((not $toprint) and $overhead);
        return undef unless @$toprint;
        
        my $entry=$self->_getEntry($toprint);
        return $entry;
    }
    
    sub _getEntry
    {
        my $self=shift;
        my $top=shift;
        die "to print is emtpy" unless @$top;
        my $chr=$top->[0]{chr};
        my $pos=$top->[0]{pos};
        my $rc=$top->[0]{refc};
        
        foreach my $t (@$top)
        {
            die "reference ids do not agree within the synchronized entry" unless $t->{chr} eq $chr;
            die "positions do not agree within the synchronized entry" unless $t->{pos} eq $pos;
            #die "reference characters do not agree within one synchronized entry" unless $t->{refc} eq $rc;
        }
        
        
        return {
          chr=>$chr,
          pos=>$pos,
          refc=>$rc,
          data=>$top
        };
    }
    
    sub _splitCurrentLines
    { # sort and take the last entries
        my $self=shift;
        my $currentLines=$self->{curlines};
        
        
        $currentLines=  [ sort {$a->{chr} cmp $b->{chr}
                                or $a->{pos} <=> $b->{pos} } @$currentLines];
        
        my $toprint=[];
        my $overhead=[];
        
        my $c=$currentLines->[0];
        for(my $i=0; $i<@$currentLines; $i++)
        {
            my $e=$currentLines->[$i];
            if($c->{chr} eq $e->{chr} && $c->{pos}==$e->{pos})
            {
                push @$toprint,$e;
            }
            else
            {
                push @$overhead,$e;
            }
        }
        return ($overhead,$toprint);
    }
    
    sub _getNewLines
    {
        my $self=shift;
        # read from the files from which lines have been printed
        my $pileups=$self->{pileupfiles};
        my $toprint=shift;
        
        my $novelLines=[];
        foreach(@$toprint)
        {
            my $index =$_->{index};
            my $line=$pileups->[$index]->nextLine();
            push @$novelLines,$line if $line;
        }
        return $novelLines;
    }
    
    
    
}







{
    package PileupFile;
    use FindBin qw($RealBin);
    use lib "$RealBin/Modules";
    use Pileup;
    my $counter=0;
    sub new
    {
        my $class=shift;
        my $path=shift;
        my $pp=shift;
        open my $fh,"<$path" or die "Could not open path";
        $counter||=0;
        
        my $e={
            file=>$path,
            fh=>$fh,
            pp=>$pp,
            index=>$counter
        };
        $counter++;
        return bless $e, __PACKAGE__;
        
    }
    
    sub nextLine
    {
        my $self=shift;
        my $fh=$self->{fh};
        my $index=$self->{index};
        my $pp=$self->{pp};
        
        my $line=<$fh>;
        return undef unless $line;
        chomp $line;

        
        my $entry= $pp->($line);
        $entry->{index}=$index;
        
        return $entry;
        
    }

}


{
    package SyncTest;
    use strict;
    use warnings;
    use FindBin qw($RealBin);
    use lib "$RealBin/Modules";
    use Test;
    use Test::PileupParser;
    use Pileup;

    
    
    sub _getPileupSliderForString
    {
        my $str=shift;
        my $minqual=shift;
        my $index=shift;
        my $pp=get_pileup_parser("illumina",1,1,2000000000,$minqual);


        open my $ofh,"<",\$str or die "could not open string filehandle";
        
        my $cr=bless {

            fh=>$ofh,
            index=>$index,
            pp=>$pp

        },"PileupFile";
        return $cr;
    }
    
    sub runTests
    {

        #runVariabilityTests();
        run_PileupParserTests();
        test_pileupreader();
        test_synchronization();
        
        exit;
    }
    
    sub test_synchronization
    {
        # 2L 2LHet 2R 2RHet 3L 3LHet 3R 3RHet 4 U_minus_mitoch Uextra X XHet YHet dmel_mitochondrion_genome gi|42519920|ref|NC_002978.6
        
        my ($s1,$s2,$s3);
        my $pfs;
        my $s;
        my $e;
        $s1=
        "2L\t10\tA\t1\tC\tT\n".
        "2L\t11\tA\t8\tA\tT\n".
        "2L\t12\tA\t5\tG\tT\n";
        $s2=
        "2L\t9\tA\t1\tC\tT\n".
        "2L\t11\tA\t8\tA\tT\n".
        "2L\t12\tA\t5\tG\tT\n";
        $s3=
        "2L\t11\tA\t1\tC\tT\n".
        "2L\t12\tA\t8\tA\tT\n".
        "2L\t13\tA\t5\tG\tT\n";
        $pfs= [_getPileupSliderForString($s1,20,0),_getPileupSliderForString($s2,20,1),_getPileupSliderForString($s3,20,2)];
        $s=Synchronizer->new($pfs);
        
        $e=$s->next();
        is($e->{chr},"2L","test synchronization; chromosome correct");
        is($e->{pos},9,"test synchronization; position correct");
        is($e->{refc},"A","test synchronization; ref-character is correct");
        is($e->{data}[0]{index},1,"index is correct");
        is(scalar(@{$e->{data}}),1,"test synchronization; number of entries is correct");
        $e=$s->next();
        is($e->{pos},10,"test synchronization; position correct");
        is(scalar(@{$e->{data}}),1,"test synchronization; number of entries is correct");
        is($e->{data}[0]{index},0,"index is correct");
        $e=$s->next();
        is($e->{pos},11,"test synchronization; position correct");
        is(scalar(@{$e->{data}}),3,"test synchronization; number of entries is correct");
        $e=$s->next();
        is($e->{pos},12,"test synchronization; position correct");
        is(scalar(@{$e->{data}}),3,"test synchronization; number of entries is correct");
        $e=$s->next();
        is($e->{pos},13,"test synchronization; position correct");
        is(scalar(@{$e->{data}}),1,"test synchronization; number of entries is correct");
        is($e->{chr},"2L","test synchronization; chromosome correct");
        is($e->{data}[0]{index},2,"index is correct");
        $e=$s->next();
        not_exists($e,"test synchronization; correct no more windows");
        
        # different chromosomes
        # 2L 2LHet 2R 2RHet 3L 3LHet 3R 3RHet 4 U_minus_mitoch Uextra X XHet YHet 
        $s1=
        "2L\t10\tA\t1\tC\tT\n".
        "YHet\t8\tA\t8\tA\tT\n".
        "YHet\t12\tA\t5\tG\tT\n";
        $s2=
        "2R\t9\tA\t1\tC\tT\n".
        "3R\t10\tA\t8\tA\tT\n".
        "3R\t12\tA\t5\tG\tT\n";
        $s3=
        "3R\t11\tA\t1\tC\tT\n".
        "YHet\t8\tA\t8\tA\tT\n".
        "YHet\t12\tA\t5\tG\tT\n";
        $pfs= [_getPileupSliderForString($s1,20,0),_getPileupSliderForString($s2,20,1),_getPileupSliderForString($s3,20,2)];
        $s=Synchronizer->new($pfs);
        $e=$s->next();
        is($e->{chr},"2L","test synchronization; chromosome correct");
        is($e->{pos},10,"test synchronization; position correct");
        is(scalar(@{$e->{data}}),1,"test synchronization; number of entries is correct");
        $e=$s->next();
        is($e->{chr},"2R","test synchronization; chromosome correct");
        is($e->{pos},9,"test synchronization; position correct");
        is(scalar(@{$e->{data}}),1,"test synchronization; number of entries is correct");
        $e=$s->next();
        is($e->{chr},"3R","test synchronization; chromosome correct");
        is($e->{pos},10,"test synchronization; position correct");
        is(scalar(@{$e->{data}}),1,"test synchronization; number of entries is correct");
        $e=$s->next();
        is($e->{chr},"3R","test synchronization; chromosome correct");
        is($e->{pos},11,"test synchronization; position correct");
        is(scalar(@{$e->{data}}),1,"test synchronization; number of entries is correct");
        $e=$s->next();
        is($e->{chr},"3R","test synchronization; chromosome correct");
        is($e->{pos},12,"test synchronization; position correct");
        is(scalar(@{$e->{data}}),1,"test synchronization; number of entries is correct");
        $e=$s->next();
        is($e->{chr},"YHet","test synchronization; chromosome correct");
        is($e->{pos},8,"test synchronization; position correct");
        is(scalar(@{$e->{data}}),2,"test synchronization; number of entries is correct");
        $e=$s->next();
        is($e->{chr},"YHet","test synchronization; chromosome correct");
        is($e->{pos},12,"test synchronization; position correct");
        is(scalar(@{$e->{data}}),2,"test synchronization; number of entries is correct");
        $e=$s->next();
        not_exists($e,"test synchronization; correct no more windows");
        
    }
    
    sub test_pileupreader
    {
        my $str;
        my $pilsl;
        my $e;
        
        $str=
        "2L\t1\tA\t9\tCCCCCAAAA\tTTTTTTTTT\n".
        "2L\t2\tA\t8\tCCCC*AAA\tTTTTTTTT\n".
        "2L\t3\tA\t5\tGGGTTNN\tTTTTTTT\n"
        ;

        $pilsl=_getPileupSliderForString($str,20,3);  # win, step, mincount, mincov, minqual
        $e=$pilsl->nextLine();
        
        is($e->{A},4,"PileupFile; A count ok");
        is($e->{index},3,"PileupFile; index is ok");
        is($e->{pos},1,"PileupFile; position is ok");
        is($e->{refc},'A',"PileupFile; reference character is ok");
        is($e->{iscov},1,"PileupFile; is covered is ok");
        is($e->{del},0,"PileupeFile; deletions are ok");
        is($e->{ispuresnp},1,"PileupFile; is pure snp is ok");
        is($e->{totcov},9,"PileupFile; total coverage is ok");
        is($e->{N},0,"PileupFile; N count is ok");
        is($e->{chr},"2L","PileupFile; contig is ok");
        
        $e=$pilsl->nextLine();
        is($e->{A},3,"PileupFile; A count ok");
        is($e->{index},3,"PileupFile; index is ok");
        is($e->{pos},2,"PileupFile; position is ok");
        is($e->{del},1,"PileupeFile; deletions are ok");
        
        $e=$pilsl->nextLine();
        is($e->{A},0,"PileupFile; A count ok");
        is($e->{index},3,"PileupFile; index is ok");
        is($e->{pos},3,"PileupFile; position is ok");
        is($e->{del},0,"PileupeFile; deletions are ok");
        is($e->{N},2,"PileupFile; N count is ok");
        

    }
    
    sub test_sorting
    {
        
    }
}








=head1 NAME

synchronize-pileup.pl - Synchronizes multiple pileup files and reports a lightwight statistic for every position (deprecated: we recommend the use of mpileup2sync.pl instead)

=head1 SYNOPSIS

 synchronize-pileup.pl --input pop1.pileup --input pop2.pileup --min-qual 20 --output output.txt

=head1 OPTIONS

=over 4

=item B<--input>

The input file in pileup format. At least two have to be specified. Mandatory parameter

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

a sorted pileup file; Sorting has to be done by (i) ascii reference-id and (ii) numerical position in the ref-seq
This can be achieved using the Unix-command:
 
 sort -k 1,1 -k 2,2n unsorted.pileup > sorted.pileup
 
 example of correct sorting:
 2L 2LHet 2R 2RHet 3L 3LHet 3R 3RHet 4 U_minus_mitoch Uextra X XHet YHet dmel_mitochondrion_genome gi|42519920|ref|NC_002978.6
 

=head2 OUTPUT

Output is a single tab delimited file which contains a lightwight representation of every pileup file.
Every pileup file represents a population and will be parsed into a list of A-count:T-count:C-count:G-count:N-count:*-count

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


=head2 CONTROL

You may check the sorting of the pileup file using the command:
 
 cat sorted | awk '{print $1}' | uniq

In the end please make sure that the synchronized file is correctly synchronized.
The corresponding pileup entries from one file should be again using the command

 cat synchronized-output | awk '{print $1}' | uniq

This command provides a list of reference sequences. The file has been correctly synchronized if every reference sequence in this list occures exactly once!

=head2 USAGE NOTE

This script was created before the mpileup file format was introduced. Currently we strongly recommend to use C<mpileup2sync.pl> instead.
However, in case only the pileup files are available this script may still be useful.

=head1 AUTHORS

Robert Kofler

Christian Schloetterer

=cut
