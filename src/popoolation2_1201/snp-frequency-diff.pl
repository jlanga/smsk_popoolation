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


my $input;
my $outputprefix="";
my $help=0;
my $test=0;

my $mincount=2;
my $mincoverage=4;
my $usermaxcoverage;
my $regionEncoded="";


GetOptions(
    "input=s"	        =>\$input,
    "output-prefix=s"   =>\$outputprefix,
    "min-count=i"       =>\$mincount,
    "min-coverage=i"    =>\$mincoverage,
    "max-coverage=s"    =>\$usermaxcoverage,
    "region=s"          =>\$regionEncoded,
    "test"              =>\$test,
    "help"	        =>\$help
) or pod2usage(-msg=>"Wrong options",-verbose=>1);

pod2usage(-verbose=>2) if $help;
SFDTest::runTests() if $test;

pod2usage(-msg=>"Input file does not exist",-verbose=>1) unless -e $input;
pod2usage(-msg=>"No output file has been provided",-verbose=>1) unless $outputprefix;
pod2usage(-msg=>"Minimum coverage <1 not allowed",-verbose=>1) if $mincoverage<1;
pod2usage(-msg=>"Provide a maximum coveage",-verbose=>1) unless $usermaxcoverage;
pod2usage(-msg=>"Minimumc count needs to be at least 1") unless $mincount;

my $paramfile=$outputprefix.".params";
open my $pfh, ">",$paramfile or die "Could not open $paramfile\n";
print $pfh "Using input\t$input\n";
print $pfh "Using outputprefix\t$outputprefix\n";
print $pfh "Using min-count\t$mincount\n";
print $pfh "Using min-coverage\t$mincoverage\n";
print $pfh "Using max-coverage\t$usermaxcoverage\n";
print $pfh "Using region\t$regionEncoded\n";
print $pfh "Using test\t$test\n";
print $pfh "Using help\t$help\n";
close $pfh;

my($output_rc,$output_pw)=($outputprefix."_rc",$outputprefix."_pwc");
open my $ifh, "<", $input or die "Could not open input file";
open my $ofhrc, ">", $output_rc or die "Could not open output file";
open my $ofhpw, ">", $output_pw or die "Could not open output file";

# parse the region if provided
my $region=undef;
$region=Utility::parse_region($regionEncoded) if $regionEncoded;


my $maxcoverage=get_max_coverage($input,$usermaxcoverage);
my $parser=Utility::get_lightwightparser($mincount,$mincoverage,$maxcoverage);

# header definition
my $headerwriten=0;
my $samples =0;

while(my $line = <$ifh>)
{
    chomp $line;
    #first write the header
    unless($headerwriten)
    {
        my @data=split /\s+/,$line;
        shift @data; shift @data; shift @data;
        $samples=@data;
        my($header_rc,$header_pw)=Utility::get_header($samples);
        print $ofhrc $header_rc."\n";
        print $ofhpw $header_pw."\n";
        $headerwriten=1;
    }
    
    # second define a region
    if($region)
    {
        my ($chr,$pos)=split /\s+/,$line;
        next unless $chr==$region->{chr};
        next if($pos<$region->{start} or $pos>$region->{end});
    }
    
    my $hit=$parser->($line);
    # chr, pos, rc, delsum, data, issnp, isrcsnp, ispopsnp
    next unless $hit->{issnp};
    
    # print
    # chr, pos, rc, snptype[rc,pop,rc|pop], delsum, T:1.000:27/29 ...  A   0.112:17/19:67/70   0.222:30/35:40/45 
    
    my($pwc,$pwcomp)=Utility::get_pairwise_comp($hit,$samples);
    $hit->{pwchar}=$pwc; # pairwise character
    $hit->{pwcomp}=$pwcomp; # pairwise comparisions
    
    my($toprint_rc,$toprint_pw)=Utility::formatOutput($hit);
    print $ofhrc $toprint_rc."\n";
    print $ofhpw $toprint_pw."\n";
}

exit;


{
    package Utility;
    use strict;
    use warnings;
    use FindBin qw/$RealBin/;
    use lib "$RealBin/Modules";
    use Synchronized;    
    
    sub formatOutput
    {
        # rc: chr, pos, rc, alleles, alstring, delsum, snp-type, consstring, subconsstring, consfreq, subconsfreq
        # rc: chr, pos, rc, alleles, alstring, delsum, snp-type, most_diverged_char, pwfreqdif
        my $hit=shift; 
        my $toprint_rc="$hit->{chr}\t$hit->{pos}\t$hit->{refchar}\t$hit->{alleles}\t$hit->{alstring}\t$hit->{delsum}\t";
        my $toprint_pw="$hit->{chr}\t$hit->{pos}\t$hit->{refchar}\t$hit->{alleles}\t$hit->{alstring}\t$hit->{delsum}\t";
        
        # chr, position, reference character
        
        
        #snp type
        if($hit->{isrcsnp} and $hit->{ispopsnp})
        {
            $toprint_rc.="rc|pop";
            $toprint_pw.="rc|pop";
        }
        elsif($hit->{isrcsnp})
        {
            $toprint_rc.="rc";
            $toprint_pw.="rc";
            
        }
        elsif($hit->{ispopsnp})
        {
            $toprint_rc.="pop";
            $toprint_pw.="pop";
        }
        else
        {
            die "impossible";
        }
        
        my $data=$hit->{samples};
        my $consstring="";
        my $subconsstring="";
        my @consfreq=();
        my @subconsfreq=();
        
        foreach my $d (@$data)
        {
            # A, T, C, G, N, del, eucov, totcov, valid_cov (totcov), a_desc
            my $a_desc=$d->{a_desc};
            my $valid_cov=$d->{iscov};
            my $cons=$a_desc->[0];
            my $subcons=$a_desc->[1];
            
            my $consc=$cons->{a};
            my $consf=$cons->{c};
            my $subconsc=$subcons->{a};
            my $subconsf=$subcons->{c};
            
            unless($valid_cov)
            {
                $consc="N"; $subconsc="N"; $consf=0; $subconsf=0;
            }
            $subconsc="N" unless $subconsf;
            $consc="N" unless $consf;
            $consf="$consf/$d->{eucov}";
            $subconsf="$subconsf/$d->{eucov}";
            
            $consstring.=$consc;
            $subconsstring.=$subconsc;
            push @consfreq,$consf;
            push @subconsfreq,$subconsf;
        }
        
        my $consfreqstring=join("\t",@consfreq);
        my $subconsfreqstring=join("\t",@subconsfreq);
        $toprint_rc.="\t$consstring\t$subconsstring\t$consfreqstring\t$subconsfreqstring";
        
        #pairwisecomparisions
        $toprint_pw.="\t$hit->{pwchar}";
        my $pwcomp=$hit->{pwcomp};
        foreach my $pw (@$pwcomp)
        {
            # p1,p2,dif
            my $d=$pw->{dif};
            $d=sprintf("%.3f",$d) unless $d eq "na";
            $toprint_pw.="\t$d";
        }
        return ($toprint_rc,$toprint_pw);
    }
    
    
    
    
    sub get_pairwise_comp
    {
        my $hit=shift;
        my $samples=shift;
        my $data=$hit->{samples};
        
        my @alleles=( {a=>"A",c=>0} , {a=>"T",c=>0} , {a=>"C",c=>0} , {a=>"G",c=>0} );
        
        for my $i(0..$samples-1)
        {
            for my $k($i+1..$samples-1)
            {
                my $e1=$data->[$i];
                my $e2=$data->[$k];
                next unless $e1->{iscov};
                next unless $e2->{iscov};
                
                foreach my $al (@alleles)
                {
                    my $f1=$e1->{$al->{a}}/$e1->{eucov};
                    my $f2=$e2->{$al->{a}}/$e2->{eucov};
                    my $div=abs($f1-$f2);
                    $al->{c}+=$div;
                }
            }
        }
        
        @alleles=sort {$b->{c}<=>$a->{c}} @alleles;
        my $pwc=$alleles[0]->{a};
        
        my @pwcomp=();
        for my $i(0..$samples-1)
        {
            for my $k($i+1..$samples-1)
            {
                my $e1=$data->[$i];
                my $e2=$data->[$k];
                if($e1->{iscov} && $e2->{iscov})
                {
                    my $f1=$e1->{$pwc} / $e1->{eucov};
                    my $f2=$e2->{$pwc} / $e2->{eucov};
                    my $div=abs($f1-$f2);
                    
                    
                    # p1,p2,dif
                    push @pwcomp,
                    {
                        p1=>$e1->{$pwc}."/".$e1->{eucov},
                        p2=>$e2->{$pwc}."/".$e2->{eucov},
                        dif=>$div,
                    };
                    
                }
                else
                {
                    push @pwcomp,
                    {
                        dif=>"na",p1=>"-",p2=>"-"
                    };
                }
            }
        }
        
        return $pwc,\@pwcomp;
    }
    
    sub get_header
    {
        my $samples=shift;
        # rc: chr, pos, rc, alleles, alstring, delsum, consstring, subconsstring, consfreq, subconsfreq
        # rc: chr, pos, rc, alleles, alstring, delsum, most_diverged_char, pwfreqdif
        my $header_rc="##chr\tpos\trc\tallele_count\tallele_states\tdeletion_sum\tsnp_type\tmajor_alleles(maa)\tminor_alleles(mia)";
        
        for my $i (1..$samples)
        {
            $header_rc.="\tmaa_$i";
        }
        for my $i (1..$samples)
        {
            $header_rc.="\tmia_$i";
        }

        my $header_pw="##chr\tpos\trc\tallele_count\tallele_states\tdeletion_sum\tsnp_type\tmost_variable_allele";
        for my $i(1..$samples)
        {
            for my $k($i+1..$samples)
            {
                $header_pw.="\tdiff:$i-$k";
            }
        }
        return ($header_rc,$header_pw);
    }
    
    sub parse_region
    {
        my $reg=shift;
        my($chr,$start,$end)=$reg=~m/^(\w+):(\d+)-(\d+)$/;
        
        return {
            chr=>$chr,
            start=>$start,
            end=>$end
        };
    }
    
    
    sub get_lightwightparser
    {
        my $mincount=shift;
        my $mincov=shift;
        my $maxcov=shift;
        my $sp=get_sumsnp_synparser($mincount,$mincov,$maxcov);
        
        return sub
        {
            my $line=shift;
            my $e=$sp->($line);
            # iscov, issnp, ispuresnps
            # chr, pos, refchar, minsampcov, samples (A T C G N del)
            $e=_update_lightwight($e,$mincount);
            return $e;
        }
    }
    
    sub _update_lightwight
    {
        my $entry=shift;
        my $mincount=shift;
        
        my $rc=$entry->{refchar};
        my $delcount=0;
        my $isrcsnp=0;
        my @alleles=({a=>"A",c=>0},{a=>"T",c=>0},{a=>"C",c=>0},{a=>"G",c=>0});
        foreach my $s (@{$entry->{samples}})
        {

            my @cons=(); push @cons,{a=>"A",c=>$s->{A} };  push @cons,{a=>"T",c=>$s->{T} };  push @cons,{a=>"C",c=>$s->{C} };  push @cons,{a=>"G",c=>$s->{G}};
            @cons=sort { $b->{c} <=> $a->{c} } @cons;
            $s->{a_desc} = \@cons;
            
            # perform the following operations only when the data have the correct coverage
            if($s->{iscov})
            {
                 my $delcount+=$s->{del};
                # refcharsnp
                $isrcsnp=1 if $s->{a_desc}[0]{a} ne $rc;
                
                # the alleles
                foreach my $al (@alleles)
                {
                    $al->{c}+=$s->{$al->{a}};
                }
            }

        }
        @alleles=sort {$b->{c}<=>$a->{c}} @alleles;
        my @valalles=grep {$_->{c}>=$mincount} @alleles;
        my $validalleles=@valalles;
        my $alstring=join("/",map {$_->{a}} @valalles);
        
        my $ispopsnp=$entry->{issnp};
        $isrcsnp=0 if $rc eq "N";
        my $issnp=($ispopsnp or $isrcsnp)? 1:0;
        $entry->{delsum}=$delcount;
        $entry->{alleles}=$validalleles;
        $entry->{alstring}=$alstring;
        $entry->{ispopsnp}=$ispopsnp;
        $entry->{isrcsnp}=$isrcsnp;
        $entry->{issnp}=$issnp;
        return $entry;
    }
    



    

}


{
    package SFDTest;
    use strict;
    use warnings;
    use FindBin qw/$RealBin/;
    use lib "$RealBin/Modules";
    use Test::TMaxCoverage;
    use Test::TSynchronized;
    use Test;
    
    sub runTests
    {
        run_MaxCoverageTests();
        run_SynchronizedTests();
        test_parseLightwight();
        test_pairwiseComparisions();
        exit;
    }
    
    sub test_pairwiseComparisions
    {
        my($pch,$p);
        my $parser=Utility::get_lightwightparser(3,4,[10,10]);
        ($pch,$p)=Utility::get_pairwise_comp($parser->("2L\t1\tC\t0:1:4:0:0:0\t0:0:4:0:0:0"),2);
        is($pch,"T","pairwise comparision: pairwise character correct");
        is($p->[0]{dif},0.2,"pairwise comparision: difference is correct");
        is($p->[0]{p1},"1/5","pairwise comparision: allele count is correct");
        is($p->[0]{p2},"0/4","pairwise comparision: allele count is correct");
        
        ($pch,$p)=Utility::get_pairwise_comp($parser->("2L\t1\tC\t0:0:4:0:0:0\t1:3:4:0:0:0"),2);
        is($pch,"C","pairwise comparision: pairwise character correct");
        is($p->[0]{dif},0.5,"pairwise comparision: difference is correct");
        is($p->[0]{p1},"4/4","pairwise comparision: allele count is correct");
        is($p->[0]{p2},"4/8","pairwise comparision: allele count is correct");
        
        $parser=Utility::get_lightwightparser(3,4,[10,10,10]);
        ($pch,$p)=Utility::get_pairwise_comp($parser->("2L\t1\tC\t0:0:4:0:0:0\t0:1:4:2:0:0\t1:3:4:0:0:0"),3);
        is($pch,"C","pairwise comparision: pairwise character correct");
        ok(abs($p->[0]{dif}-0.42857)<0.001,"pairwise comparision: difference is correct");
        is($p->[0]{p1},"4/4","pairwise comparision: allele count is correct");
        is($p->[0]{p2},"4/7","pairwise comparision: allele count is correct");
        ok(abs($p->[1]{dif}-0.5)<0.001,"pairwise comparision: difference is correct");
        is($p->[1]{p1},"4/4","pairwise comparision: allele count is correct");
        is($p->[1]{p2},"4/8","pairwise comparision: allele count is correct");
        ok(abs($p->[2]{dif}-0.0714)<0.001,"pairwise comparision: difference is correct");
        is($p->[2]{p1},"4/7","pairwise comparision: allele count is correct");
        is($p->[2]{p2},"4/8","pairwise comparision: allele count is correct");
        
    }
    
    sub test_parseLightwight
    {
        my $r;
        my $parser=Utility::get_lightwightparser(3,5,[10,10]);
        # $line, $mincount, $mincov, $maxcov
        $r=$parser->("2L\t1\tC\t0:1:4:0:0:0\t0:2:3:0:0:0");
        is($r->{issnp},1,"parseLightwight: is a SNP correct");
        is($r->{ispopsnp},1,"parseLightwight: is a population SNP correct");
        is($r->{isrcsnp},0,"parseLightwight: is a SNP with the reference character: correct");
        is($r->{chr},"2L","parseLightwight: chromosome correct");
        is($r->{refchar},"C","parseLightwight: reference character is correct");
        is($r->{pos},1,"parseLightwight: position is correct");
        is($r->{alleles},2,"parseLightwight; Number of alleles correct");
        is($r->{alstring},"C/T","parseLightwight; Allelestring correct");
        
        # an N in the reference should not fuck up the result
        $r=$parser->("2L\t1\tN\t0:1:4:0:0:0\t0:2:3:0:0:0");
        is($r->{issnp},1,"parseLightwight: is a SNP correct");
        is($r->{ispopsnp},1,"parseLightwight: is a population SNP correct");
        is($r->{isrcsnp},0,"parseLightwight: is a SNP with the reference character: correct");
        is($r->{chr},"2L","parseLightwight: chromosome correct");
        is($r->{refchar},"N","parseLightwight: reference character is correct");
        is($r->{pos},1,"parseLightwight: position is correct");
        is($r->{alleles},2,"parseLightwight; Number of alleles correct");
        is($r->{alstring},"C/T","parseLightwight; Allelestring correct");
        
        # coverage in one population to low?
        $r=$parser->("2L\t1\tN\t0:0:4:0:0:0\t0:2:3:0:0:0");
        is($r->{issnp},0,"parseLightwight: is a SNP correct");
        is($r->{ispopsnp},0,"parseLightwight: is a population SNP correct");
        is($r->{isrcsnp},0,"parseLightwight: is a SNP with the reference character: correct");
        is($r->{chr},"2L","parseLightwight: chromosome correct");
        is($r->{refchar},"N","parseLightwight: reference character is correct");
        is($r->{pos},1,"parseLightwight: position is correct");
        is($r->{alleles},1,"parseLightwight; Number of alleles correct");
        is($r->{alstring},"C","parseLightwight; Allelestring correct");
        
        
        # NO SNP
        $r=$parser->("2L\t2\tC\t0:0:5:0:0:0\t0:0:5:0:0:0");
        is($r->{issnp},0,"parseLightwight: is a SNP correct");
        is($r->{ispopsnp},0,"parseLightwight: is a population SNP correct");
        is($r->{isrcsnp},0,"parseLightwight: is a SNP with the reference character: correct");
        is($r->{chr},"2L","parseLightwight: chromosome correct");
        is($r->{refchar},"C","parseLightwight: reference character is correct");
        is($r->{pos},2,"parseLightwight: position is correct");
        is($r->{alleles},1,"parseLightwight; Number of alleles correct");
        is($r->{alstring},"C","parseLightwight; Allelestring correct");
        
        
        # a reference character SNP
        $r=$parser->("3R\t2\tC\t0:5:0:0:0:0\t0:5:0:0:0:0");
        is($r->{issnp},1,"parseLightwight: is a SNP correct");
        is($r->{ispopsnp},0,"parseLightwight: is a population SNP correct");
        is($r->{isrcsnp},1,"parseLightwight: is a SNP with the reference character: correct");
        is($r->{chr},"3R","parseLightwight: chromosome correct");
        is($r->{refchar},"C","parseLightwight: reference character is correct");
        is($r->{pos},2,"parseLightwight: position is correct");
        is($r->{alleles},1,"parseLightwight; Number of alleles correct");
        is($r->{alstring},"T","parseLightwight; Allelestring correct");
        
        $r=$parser->("3R\t2\tC\t0:5:0:0:0:0\t0:5:2:2:2:2");
        is($r->{issnp},1,"parseLightwight: is a SNP correct");
        is($r->{ispopsnp},0,"parseLightwight: is a population SNP correct");
        is($r->{isrcsnp},1,"parseLightwight: is a SNP with the reference character: correct");
        is($r->{chr},"3R","parseLightwight: chromosome correct");
        is($r->{refchar},"C","parseLightwight: reference character is correct");
        is($r->{pos},2,"parseLightwight: position is correct");
        is($r->{alleles},1,"parseLightwight: Number of alleles correct");
        is($r->{alstring},"T","parseLightwight: Allelestring correct");
        is($r->{samples}[1]{eucov},9,"parseLightwight: eu-coverage is correct");
        is($r->{samples}[1]{totcov},13,"parseLightwight: tot-coverage is correct");
        is($r->{samples}[1]{del},2,"parseLightwight: del count is correct");    
        is($r->{samples}[1]{N},2,"parseLightwight: N count is correct");
        is($r->{samples}[1]{A},0,"parseLightwight: A count is correct");
        is($r->{samples}[1]{G},2,"parseLightwight: G count is correct");
        is($r->{samples}[1]{C},2,"parseLightwight: C count is correct");
        is($r->{samples}[1]{T},5,"parseLightwight: T count is correct");
        
        # counts
        $r=$parser->("3R\t2\tC\t0:5:0:0:0:0\t5:4:3:2:1:6");
        is($r->{samples}[1]{eucov},14,"parseLightwight: eu-coverage is correct");
        is($r->{samples}[1]{totcov},21,"parseLightwight: tot-coverage is correct");
        is($r->{samples}[1]{del},6,"parseLightwight: del count is correct");    
        is($r->{samples}[1]{N},1,"parseLightwight: N count is correct");
        is($r->{samples}[1]{A},5,"parseLightwight: A count is correct");
        is($r->{samples}[1]{G},2,"parseLightwight: G count is correct");
        is($r->{samples}[1]{C},3,"parseLightwight: C count is correct");
        is($r->{samples}[1]{T},4,"parseLightwight: T count is correct");
        
        
        # reference character is a 'N'
        $r=$parser->("3R\t2\tN\t0:5:0:0:0:0\t0:5:0:0:0:0");
        is($r->{issnp},0,"parseLightwight: is a SNP correct");
        is($r->{ispopsnp},0,"parseLightwight: is a population SNP correct");
        is($r->{isrcsnp},0,"parseLightwight: is a SNP with the reference character: correct");
        is($r->{chr},"3R","parseLightwight: chromosome correct");
        is($r->{refchar},"N","parseLightwight: reference character is correct");
        is($r->{pos},2,"parseLightwight: position is correct");
        is($r->{alleles},1,"parseLightwight; Number of alleles correct");
        is($r->{alstring},"T","parseLightwight; Allelestring correct");
        
        
        # a pop and rc SNP
        $r=$parser->("3R\t2\tC\t0:5:4:0:0:0\t3:5:0:0:0:0");
        is($r->{issnp},1,"parseLightwight: is a SNP correct");
        is($r->{ispopsnp},1,"parseLightwight: is a population SNP correct");
        is($r->{isrcsnp},1,"parseLightwight: is a SNP with the reference character: correct");
        is($r->{chr},"3R","parseLightwight: chromosome correct");
        is($r->{refchar},"C","parseLightwight: reference character is correct");
        is($r->{pos},2,"parseLightwight: position is correct");
        is($r->{alleles},3,"parseLightwight; Number of alleles correct");
        is($r->{alstring},"T/C/A","parseLightwight; Allelestring correct");
        
    }
}


=head1 NAME

SNP-frequency-diff.pl - Identify differences in allele frequencies between at least two populations

=head1 SYNOPSIS

 SNP-frequency-diff.pl --input populations.sync --output-prefix pop_diff --min-count 2 --min-coverage 4  --max-coverage 2%

=head1 OPTIONS

=over 4

=item B<--input>

The input file. Has to be synchronized pileup file. Mandatory parameter

=item B<--output-prefix>

The prefix of the output files. Mandatory parameter

The script will create two files:

 "output_prefix"_pwc: contains the differences in allele frequencies for all pairwise comparisions
 "output_prefix"_rc: shows the major and minor allele in a succinct format for all populations

=item B<--region>

provide this option if only a certain subregion of a contig should be analysed; has to be of the form "contig:start-end" for example: "chr2:10000-20000"; default=empty

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

=item B<--test>

Run the unit tests for this script. 

=item B<--help>

Display help for this script

=back

=head1 DETAILS


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

=head2 OUTPUT <output-prefix>_rc

 ##chr   pos     rc      allele_count    allele_states   deletion_sum    snp_type        major_alleles(maa)      minor_alleles(mia)      maa_1   maa_2   maa_3   maa_4   maa_5   mia_1   mia_2   mia_3   mia_4   mia_5
 gi|42519920|ref|NC_002978.6     3530    A       2       A/G     0       pop     AAAAA   GGNNG   54/55   76/78   35/35   27/27   17/18   1/55    2/78    0/35    0/27    1/18
 gi|42519920|ref|NC_002978.6     5450    T       2       T/A     0       pop     TTTTT   ANANN   51/54   96/96   9/10    33/33   23/23   3/54    0/96    1/10    0/33    0/23


 col1: reference chromosome (contig)
 col2: referenced position
 col3: reference character
 col4: number of alleles found in all populations
 col5: allele characters in all populations (sorted by counts in all populations)
 col6: sum of deletions in all populations (should be zero, if not the postion may not be reliable)
 col7: SNP type: [pop, rc, rc|pop]; pop.. a SNP within or between the populations; rc.. a SNP between the reference sequence character and the consensus of at least one populaton; rc|pop..both
 col8: most frequent allele in all populations [12345..]
 col9: second most frequent allele in all populations [12345..]
 col10 - col9+n: frequencies of the most frequent allele (major) in the form "allele-count/coverage"
 col10+n - col9+2n: frequencies of the second most frequent allele (minor) in the form "allele-count/coverage"

=head2 OUTPUT <output-prefix>_pwc

 ##chr   pos     rc      allele_count    allele_states   deletion_sum    snp_type        most_variable_allele    diff:1-2        diff:1-3        diff:1-4        diff:1-5        diff:2-3        dif
 f:2-4        diff:2-5        diff:3-4        diff:3-5        diff:4-5
 gi|42519920|ref|NC_002978.6     3530    A       2       A/G     0       pop     A       0.007   0.018   0.018   0.037   0.026   0.026   0.030   0.000   0.056   0.056
 gi|42519920|ref|NC_002978.6     5450    T       2       T/A     0       pop     A       0.056   0.044   0.056   0.056   0.100   0.000   0.000   0.100   0.100   0.000
 
 col1: reference chromosome (contig)
 col2: referenced position
 col3: reference character
 col4: number of alleles found in all populations
 col5: allele characters in all populations (sorted by counts in all populations)
 col6: sum of deletions in all populations (should be zero, if not the postion may not be reliable)
 col7: SNP type: [pop, rc, rc|pop]; pop.. a SNP within or between the populations; rc.. a SNP between the reference sequence character and the consensus of at least one populaton; rc|pop..both
 col8: allele showing the most variability in terms of frequency between all pairwise comparisions; for biallelic SNPs one allele is picked at random
 diff:1-2 frequency difference for allele <col8> between population-1 and population-2
 diff:2-4 frequency difference for allele <col8> between population-2 and population-4
 ...


=head1 AUTHORS

Robert Kofler

Christian Schloetterer

=cut



