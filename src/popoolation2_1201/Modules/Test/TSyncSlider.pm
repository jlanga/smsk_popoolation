{
    package Test::TSyncSlider;
    use strict;
    use warnings;
    use FindBin qw/$Bin/;
    use lib "$Bin/..";
    use Test;
    use SyncSlider;
    use Synchronized;
    require Exporter;
    our @ISA = qw(Exporter);
    our @EXPORT=qw(run_SyncSliderTests);
    our @EXPORT_OK = qw();
    
    
    sub run_SyncSliderTests
    {
        testSyncSlider();
    }
    
    
    sub _getSyncSliderForString #str, window, step, mincount, mincov, maxcov
    {
        my $str=shift;
        my $window=shift;
        my $step=shift;
        my $mincount=shift;
        my $mincov=shift;
        my $maxcov=shift;
        die "No maximum coverage" unless $maxcov;
        my $sp=get_sumsnp_synparser($mincount,$mincov,$maxcov);

        open my $ofh,"<",\$str or die "could not open string filehandle";
        my $cr=bless {
            lower=>0,
            upper=>$window,
            window=>$window,
            step=>$step,
            fh=>$ofh,
            sp=>$sp,
            curwin=>[],
            buffer=>[]
        },"SyncSlider";
        return $cr;
    }
    
    sub testSyncSlider
    {
        my $teststr=
        "chr1\t1\tN\t0:0:6:0:0:0\t0:0:0:7:0:0\t0:0:8:0:0:0\n".
        "chr1\t2\tN\t0:4:8:0:0:0\t0:0:8:0:0:0\t0:0:8:0:0:0\n".
        "chr1\t3\tN\t0:0:3:0:0:0\t0:0:0:8:0:0\t0:0:13:0:0:0\n".
        "chr1\t4\tN\t0:0:5:0:0:0\t0:0:4:0:0:0\t0:0:6:0:0:0\n".
        "chr1\t5\tN\t0:0:8:0:0:0\t0:0:0:8:0:0\t0:0:8:0:0:0\n".
        "chr2\t3\tN\t0:0:11:0:0:0\t0:0:0:8:0:0\t0:0:20:0:0:0\n".
        "chr3\t1\tN\t0:0:8:0:0:0\t0:0:8:0:0:0\t0:0:8:0:0:0\n".
        "chr3\t2\tN\t0:0:1:0:0:0\t0:0:0:5:0:0\t0:0:8:0:0:0\n".
        "chr4\t4\tN\t0:0:8:0:0:0\t0:0:8:0:0:0\t0:0:8:0:0:0\n";
        my $bpsl=_getSyncSliderForString($teststr,3,1,2,4,[10000,10000,10000]);
        
        my $w=$bpsl->nextWindow();
        
        is($w->{chr},"chr1","test Synclider, correct chromosome");
        is($w->{count_covered},2,"test SyncSlider, correct number of sufficiently covered regions");
        is($w->{countpuresnp},2,"test SyncSlider, correct number of snps in region");
        is($w->{start},0,"test SyncSlider, correct start position");
        is($w->{end},3,"test SyncSlider, correct end position");
        is($w->{window},3,"test SyncSlider, correct window length");
        is(scalar(@{$w->{data}}),3,"Correct number of data entries");
        
        
        
        $w=$bpsl->nextWindow();
        is($w->{chr},"chr1","test SyncSlider, correct chromosome");
        is($w->{count_covered},2,"test SyncSlider, correct number of sufficiently covered regions");
        is($w->{countpuresnp},1,"test SyncSlider, correct number of snps in region");
        is($w->{start},1,"test SyncSlider, correct start position");
        is($w->{end},4,"test SyncSlider, correct end position");
        is($w->{middle},3,"test SyncSlider, correct end position");
        is($w->{window},3,"test SyncSlider, correct window length");
        is(scalar(@{$w->{data}}),3,"Correct number of data entries");
        
        
        $w=$bpsl->nextWindow();
        is($w->{chr},"chr1","test SyncSlider, correct chromosome");
        is($w->{count_covered},2,"test SyncSlider, correct number of sufficiently covered regions");
        is($w->{countpuresnp},1,"test SyncSlider, correct number of snps in region");
        is($w->{start},2,"test SyncSlider, correct start position");
        is($w->{end},5,"test SyncSlider, correct end position");
        is($w->{middle},4,"test SyncSlider, correct end position");
        is($w->{window},3,"test SyncSlider, correct window length");
        is(scalar(@{$w->{data}}),3,"Correct number of data entries");

        
        # switch to new chromosome
        $w=$bpsl->nextWindow();
        is($w->{chr},"chr2","test SyncSlider, correct chromosome");
        is($w->{count_covered},1,"test SyncSlider, correct number of sufficiently covered regions");
        is($w->{countpuresnp},1,"test SyncSlider, correct number of snps in region");
        is($w->{start},0,"test SyncSlider, correct start position");
        is($w->{end},3,"test SyncSlider, correct end position");
        is($w->{middle},2,"test SyncSlider, correct end position");
        is($w->{window},3,"test SyncSlider, correct window length");
        
        is(scalar(@{$w->{data}}),1,"Correct number of data entries");

        
        $w=$bpsl->nextWindow();
        is($w->{chr},"chr3","test SyncSlider, correct chromosome");
        is($w->{count_covered},1,"test SyncSlider, correct number of sufficiently covered regions");
        is($w->{countpuresnp},0,"test SyncSlider, correct number of snps in region");
        is($w->{start},0,"test SyncSlider, correct start position");
        is($w->{end},3,"test SyncSlider, correct end position");
        is($w->{middle},2,"test SyncSlider, correct end position");
        is($w->{window},3,"test SyncSlider, correct window length");
        is(scalar(@{$w->{data}}),2,"Correct number of data entries");

        
        $w=$bpsl->nextWindow();
        is($w->{chr},"chr4","test SyncSlider, correct chromosome");
        is($w->{count_covered},0,"test SyncSlider, correct number of sufficiently covered regions");
        is($w->{countpuresnp},0,"test SyncSlider, correct number of snps in region");
        is($w->{start},0,"test SyncSlider, correct start position");
        is($w->{end},3,"test SyncSlider, correct end position");
        is($w->{middle},2,"test SyncSlider, correct end position");
        is($w->{window},3,"test SyncSlider, correct window length");
        is(scalar(@{$w->{data}}),0,"Correct number of data entries");
        
        $w=$bpsl->nextWindow();
        is($w->{chr},"chr4","test SyncSlider, correct chromosome");
        is($w->{count_covered},1,"test SyncSlider, correct number of sufficiently covered regions");
        is($w->{countpuresnp},0,"test SyncSlider, correct number of snps in region");
        is($w->{start},1,"test SyncSlider, correct start position");
        is($w->{end},4,"test SyncSlider, correct end position");
        is($w->{middle},3,"test SyncSlider, correct end position");
        is($w->{window},3,"test SyncSlider, correct window length");
        is(scalar(@{$w->{data}}),1,"Correct number of data entries");
        
        $w=$bpsl->nextWindow();
        not_exists($w,"correct end of file reached");
        $w=$bpsl->nextWindow();
        not_exists($w,"correct end of file reached");
        
        
        # weired characters
        $teststr=
        "chr1\t1\tN\t0:0:6:0:1:0\t0:0:0:7:3:0\t0:0:8:0:2:0\n".
        "chr1\t2\tN\t0:4:8:0:0:0\t0:0:8:0:0:0\t0:0:8:0:0:2\n".
        "chr1\t3\tN\t0:0:3:0:0:0\t0:0:0:8:0:0\t0:0:13:0:2:2\n";
        $bpsl=_getSyncSliderForString($teststr,3,1,2,4,[10000,10000,10000]);
        #str, window, step, mincount, mincov, maxcov
        
        $w=$bpsl->nextWindow();
        is($w->{chr},"chr1","test SyncSlider, correct chromosome");
        is($w->{count_covered},1,"test SyncSlider, correct number of sufficiently covered regions");
        is($w->{countpuresnp},1,"test SyncSlider, correct number of snps in region");
        is($w->{start},0,"test SyncSlider, correct start position");
        is($w->{end},3,"test SyncSlider, correct end position");
        is($w->{window},3,"test SyncSlider, correct window length");
          
  }
    
}
1;