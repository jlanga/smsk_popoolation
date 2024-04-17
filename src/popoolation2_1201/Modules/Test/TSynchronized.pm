{
    package Test::TSynchronized;
    use strict;
    use warnings;
    use FindBin qw/$Bin/;
    use lib "$Bin/..";
    use Test;
    use Synchronized;
    require Exporter;
    our @ISA = qw(Exporter);
    our @EXPORT=qw(run_SynchronizedTests);
    our @EXPORT_OK = qw();
    
    
    sub run_SynchronizedTests
    {
        test_basic_syncparser();
        test_sumsnp_syncparser();
    }
    
    
    sub test_basic_syncparser
    {
        my $p;
        my $r;
        
        $p=get_basic_syncparser();
        $r=$p->("chr1\t1\tN\t1:2:3:4:5:6\t7:8:9:10:11:12\t13:14:15:16:17:18");
        is($r->{chr},"chr1", "test basic sync parser: Chromosome is OK");
        is($r->{pos},"1","test basic sync parser: Position is OK");
        is($r->{refchar},"N","test basic sync parser: Reference character is OK");
        is(scalar(@{$r->{samples}}),3,"test basic sync parser: number of samples is OK");
        is($r->{samples}[0]{A},"1","test basic sync parser: number of A's is OK");
        is($r->{samples}[0]{T},"2","test basic sync parser: number of T's is OK");      
        is($r->{samples}[0]{C},"3","test basic sync parser: number of C's is OK");
        is($r->{samples}[0]{G},"4","test basic sync parser: number of G's is OK");
        is($r->{samples}[0]{N},"5","test basic sync parser: number of N's is OK");
        is($r->{samples}[0]{del},"6","test basic sync parser: number of deletions is OK");
        is($r->{samples}[0]{eucov},10,"test basic sync parser: coverage is OK");
        is($r->{samples}[0]{totcov},21,"test basic sync parser: total coverage is OK");
        is($r->{samples}[1]{A},"7","test basic sync parser: number of A's is OK");
        is($r->{samples}[1]{T},"8","test basic sync parser: number of T's is OK");      
        is($r->{samples}[1]{C},"9","test basic sync parser: number of C's is OK");
        is($r->{samples}[1]{G},"10","test basic sync parser: number of G's is OK");
        is($r->{samples}[1]{N},"11","test basic sync parser: number of N's is OK");
        is($r->{samples}[1]{del},"12","test basic sync parser: number of deletions is OK");
        is($r->{samples}[2]{A},"13","test basic sync parser: number of A's is OK");
        is($r->{samples}[2]{T},"14","test basic sync parser: number of T's is OK");      
        is($r->{samples}[2]{C},"15","test basic sync parser: number of C's is OK");
        is($r->{samples}[2]{G},"16","test basic sync parser: number of G's is OK");
        is($r->{samples}[2]{N},"17","test basic sync parser: number of N's is OK");
        is($r->{samples}[2]{del},"18","test basic sync parser: number of deletions is OK");
        is($r->{samples}[2]{eucov},58,"test basic sync parser: coverage is OK");
        is($r->{samples}[2]{totcov},93,"test basic sync parser: total coverage is OK");
        
        $r=$p->("robert\t12345\tA\t6:5:4:3:2:1\t-\t0:0:0:0:0:0\t");
        is($r->{chr},"robert", "test basic sync parser: Chromosome is OK");
        is($r->{pos},"12345","test basic sync parser: Position is OK");
        is($r->{refchar},"A","test basic sync parser: Reference character is OK");
        is(scalar(@{$r->{samples}}),3,"test basic sync parser: number of samples is OK");
        is($r->{samples}[0]{A},"6","test basic sync parser: number of A's is OK");
        is($r->{samples}[0]{T},"5","test basic sync parser: number of T's is OK");      
        is($r->{samples}[0]{C},"4","test basic sync parser: number of C's is OK");
        is($r->{samples}[0]{G},"3","test basic sync parser: number of G's is OK");
        is($r->{samples}[0]{N},"2","test basic sync parser: number of N's is OK");
        is($r->{samples}[0]{del},"1","test basic sync parser: number of deletions is OK");
        is($r->{samples}[1]{A},"0","test basic sync parser: number of A's is OK");
        is($r->{samples}[1]{T},"0","test basic sync parser: number of T's is OK");      
        is($r->{samples}[1]{C},"0","test basic sync parser: number of C's is OK");
        is($r->{samples}[1]{G},"0","test basic sync parser: number of G's is OK");
        is($r->{samples}[1]{N},"0","test basic sync parser: number of N's is OK");
        is($r->{samples}[1]{del},"0","test basic sync parser: number of deletions is OK");
        is($r->{samples}[2]{A},"0","test basic sync parser: number of A's is OK");
        is($r->{samples}[2]{T},"0","test basic sync parser: number of T's is OK");      
        is($r->{samples}[2]{C},"0","test basic sync parser: number of C's is OK");
        is($r->{samples}[2]{G},"0","test basic sync parser: number of G's is OK");
        is($r->{samples}[2]{N},"0","test basic sync parser: number of N's is OK");
        is($r->{samples}[2]{del},"0","test basic sync parser: number of deletions is OK");
        
        
        $r=$p->("testchr\t12\tC\t6:5:4:3:2:1\t-\t0:0:0:0:0:0\t-\t-\t-\n");
        is($r->{chr},"testchr", "test basic sync parser: Chromosome is OK");
        is($r->{pos},"12","test basic sync parser: Position is OK");
        is($r->{refchar},"C","test basic sync parser: Reference character is OK");
        is(scalar(@{$r->{samples}}),6,"test basic sync parser: number of samples is OK");
    }
    
    
    sub test_sumsnp_syncparser
    {
        my $p;
        my $r;
        
        $p=get_sumsnp_synparser(3,3,[20,20]);
        $r=$p->("chr1\t1\tN\t2:2:0:0:0:0\t2:2:0:0:0:0");
        is($r->{chr},"chr1", "test sumsnp sync parser: Chromosome is OK");
        is($r->{pos},"1","test sumsnp sync parser: Position is OK");
        is($r->{refchar},"N","test sumsnp sync parser: Reference character is OK");
        is(scalar(@{$r->{samples}}),2,"test sumsnp sync parser: number of samples is OK");
        is($r->{samples}[0]{A},"2","test sumsnp sync parser: number of A's is OK");
        is($r->{samples}[0]{T},"2","test sumsnp sync parser: number of T's is OK");      
        is($r->{samples}[0]{C},"0","test sumsnp sync parser: number of C's is OK");
        is($r->{samples}[0]{G},"0","test sumsnp sync parser: number of G's is OK");
        is($r->{samples}[0]{N},"0","test sumsnp sync parser: number of N's is OK");
        is($r->{samples}[0]{del},"0","test sumsnp sync parser: number of deletions is OK");
        is($r->{samples}[0]{eucov},4,"test sumsnp sync parser: coverage is OK");
        is($r->{samples}[0]{totcov},4,"test sumsnp sync parser: total coverage is OK");
        is($r->{issnp},1,"test sumsnp sync parser: SNP calling is OK");
        is($r->{iscov},1,"test sumsnp sync parser: coverage is OK");
        is($r->{ispuresnp},1,"test sumsnp sync parser: SNP is not tainted (by an indel)");
        
        $r=$p->("chr1\t1\tN\t1:1:0:0:0:0\t2:2:0:0:0:0");
        is($r->{issnp},0,"test sumsnp sync parser: SNP calling is OK");
        is($r->{iscov},0,"test sumsnp sync parser: coverage is OK");
        is($r->{ispuresnp},0,"test sumsnp sync parser: SNP is not tainted (by an indel)");
        $r=$p->("chr1\t1\tN\t21:0:0:0:0:0\t0:21:0:0:0:0");
        is($r->{issnp},0,"test sumsnp sync parser: SNP calling is OK");
        is($r->{iscov},0,"test sumsnp sync parser: coverage is OK");
        is($r->{ispuresnp},0,"test sumsnp sync parser: SNP is not tainted (by an indel)");
        $r=$p->("chr1\t1\tN\t3:0:0:0:0:0\t0:3:0:0:0:0");
        is($r->{issnp},1,"test sumsnp sync parser: SNP calling is OK");
        is($r->{iscov},1,"test sumsnp sync parser: coverage is OK");
        is($r->{ispuresnp},1,"test sumsnp sync parser: SNP is not tainted (by an indel)");
        
        $p=get_sumsnp_synparser(3,3,[10,20]);
        $r=$p->("chr1\t1\tN\t3:0:0:0:0:0\t0:3:0:0:0:0");
        is($r->{issnp},1,"test sumsnp sync parser: SNP calling is OK");
        is($r->{iscov},1,"test sumsnp sync parser: coverage is OK");
        is($r->{ispuresnp},1,"test sumsnp sync parser: SNP is not tainted (by an indel)");
        $p=get_sumsnp_synparser(3,3,[10,20]);
        $r=$p->("chr1\t1\tN\t11:0:0:0:0:0\t0:5:0:0:0:0");
        is($r->{issnp},0,"test sumsnp sync parser: SNP calling is OK");
        is($r->{iscov},0,"test sumsnp sync parser: coverage is OK");
        is($r->{ispuresnp},0,"test sumsnp sync parser: SNP is not tainted (by an indel)");
        $p=get_sumsnp_synparser(3,3,[20,10]);
        $r=$p->("chr1\t1\tN\t11:0:0:0:0:0\t0:5:0:0:0:0");
        is($r->{issnp},1,"test sumsnp sync parser: SNP calling is OK");
        is($r->{iscov},1,"test sumsnp sync parser: coverage is OK");
        is($r->{ispuresnp},1,"test sumsnp sync parser: SNP is not tainted (by an indel)");
        is($r->{ignore},0,"test sumsnp sync parser: ignore flag set correctly");
        
        
        # indels
        $p=get_sumsnp_synparser(3,3,[20,10]);
        $r=$p->("chr1\t1\tN\t11:0:0:0:0:3\t0:5:0:0:0:0");
        is($r->{issnp},0,"test sumsnp sync parser: SNP calling is OK");
        is($r->{iscov},0,"test sumsnp sync parser: coverage is OK");
        is($r->{ispuresnp},0,"test sumsnp sync parser: SNP is not tainted (by an indel)");
        is($r->{ignore},1,"test sumsnp sync parser: ignore flag set correctly");
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
        my $bpsl=_getSyncSliderForString($teststr,3,1,2,4);
        
        my $w=$bpsl->nextWindow();
        
        is($w->{chr},"chr1","test Synclider, correct chromosome");
        is($w->{count_covered},2,"test SyncSlider, correct number of sufficiently covered regions");
        is($w->{countpuresnp},2,"test SyncSlider, correct number of snps in region");
        is($w->{start},0,"test SyncSlider, correct start position");
        is($w->{end},3,"test SyncSlider, correct end position");
        is($w->{window},3,"test SyncSlider, correct window length");
        is(scalar(@{$w->{data}}),3,"Correct number of data entries");
        is($w->{data}[0]{ispuresnp},1,"test SyncSlider, correct identification of a pure SNP at the first position");
        is($w->{data}[1]{ispuresnp},1,"test SyncSlider, correct identification of a pure SNP at the second position");
        is($w->{data}[2]{ispuresnp},0,"test SyncSlider, correct, there is no pure SNP at position two");
        is($w->{data}[0]{issnp},1,"test SyncSlider, correct identification of a SNP at the first position");
        is($w->{data}[1]{issnp},1,"test SyncSlider, correct identification of a SNP at the second position");
        is($w->{data}[2]{issnp},0,"test SyncSlider, correct, there is no SNP at position two");
        is($w->{data}[0]{iscov},1,"test SyncSlider, first position is sufficiently covered");
        is($w->{data}[1]{iscov},1,"test SyncSlider, second position is sufficiently covered for a SNP");
        is($w->{data}[2]{iscov},0,"test SyncSlider, third position is not sufficiently covered for a SNP -> no SNP possible at this position");
        is($w->{data}[0]{pos},1,"test SyncSlider, position is correct");
        is($w->{data}[1]{pos},2,"test SyncSlider, position is correct");
        is($w->{data}[2]{pos},3,"test SyncSlider, position is correct");
        
        is($w->{data}[0]{samples}[0]{eucov},6,"test SyncSlider, coverage is correct");
        is($w->{data}[0]{samples}[1]{eucov},7,"test SyncSlider, coverage is correct");
        is($w->{data}[0]{samples}[2]{eucov},8,"test SyncSlider, coverage is correct");
        is($w->{data}[0]{samples}[0]{totcov},6,"test SyncSlider, coverage is correct");
        is($w->{data}[0]{samples}[1]{totcov},7,"test SyncSlider, coverage is correct");
        is($w->{data}[0]{samples}[2]{totcov},8,"test SyncSlider, coverage is correct");
        is($w->{data}[0]{samples}[0]{index},0,"test SyncSlider, index is correct");
        is($w->{data}[0]{samples}[1]{index},1,"test SyncSlider, index is correct");
        is($w->{data}[0]{samples}[2]{index},2,"test SyncSlider, index is correct");

        is($w->{data}[0]{samples}[0]{A},0,"test SyncSlider, count of A is correct");
        is($w->{data}[0]{samples}[0]{T},0,"test SyncSlider, count of T is correct");
        is($w->{data}[0]{samples}[0]{C},6,"test SyncSlider, count of C is correct");
        is($w->{data}[0]{samples}[0]{G},0,"test SyncSlider, count of G is correct");
        
        is($w->{data}[1]{samples}[0]{eucov},12,"test SyncSlider, coverage is ok");
        is($w->{data}[1]{samples}[0]{A},0,"test SyncSlider, count of A is correct");
        is($w->{data}[1]{samples}[0]{T},4,"test SyncSlider, count of T is correct");
        is($w->{data}[1]{samples}[0]{C},8,"test SyncSlider, count of C is correct");
        is($w->{data}[1]{samples}[0]{G},0,"test SyncSlider, count of G is correct");
        
        is($w->{data}[2]{samples}[2]{eucov},13,"test SyncSlider, coverage is ok");
        is($w->{data}[2]{samples}[2]{A},0,"test SyncSlider, count of A is correct");
        is($w->{data}[2]{samples}[2]{T},0,"test SyncSlider, count of T is correct");
        is($w->{data}[2]{samples}[2]{C},13,"test SyncSlider, count of C is correct");
        is($w->{data}[2]{samples}[2]{G},0,"test SyncSlider, count of G is correct");
        
        

    }
    
}
1;