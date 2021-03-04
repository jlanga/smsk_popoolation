{
    package Test::TMaxCoverage;
    use strict;
    use warnings;
    use FindBin qw/$Bin/;
    use lib "$Bin/..";
    use Test;
    use MaxCoverage;
    require Exporter;
    our @ISA = qw(Exporter);
    our @EXPORT  =qw(run_MaxCoverageTests);
    
    
    sub run_MaxCoverageTests
    {
        test_maxcoverage();
    }
    
    
    sub test_maxcoverage
    {
        my $file;
        my $mc;
        
        $file="chr1\t1\tN\t0:0:6:0:0:0\t0:0:0:7:0:0\t0:0:8:0:0:0\n";
        $mc=get_max_coverage(\$file,"500");
        is($mc->[0],500,"test maximum coverage method; maximum coverage is ok");
        is($mc->[1],500,"test maximum coverage script; maximum coverage is ok");
        is($mc->[2],500,"test maximum coverage script; maximum coverage is ok");
        
        $file="chr1\t1\tN\t0:0:6:0:0:0\t0:0:0:7:0:0\t0:0:8:0:0:0\n";
        $mc=get_max_coverage(\$file,"300,400,500");
        is($mc->[0],300,"test maximum coverage method; maximum coverage is ok");
        is($mc->[1],400,"test maximum coverage script; maximum coverage is ok");
        is($mc->[2],500,"test maximum coverage script; maximum coverage is ok");
        
        $file=
        "chr1\t1\tN\t1:0:0:0:0:0\t0:1:0:0:0:0\t1:0:0:0:0:0\n".
        "chr1\t2\tN\t2:0:0:0:0:0\t0:1:0:0:0:0\t1:0:0:0:0:0\n".
        "chr1\t3\tN\t3:0:0:0:0:0\t0:1:0:0:0:0\t1:0:0:0:0:0\n".
        "chr1\t4\tN\t4:0:0:0:0:0\t0:1:0:0:0:0\t1:0:0:0:0:0\n".
        "chr1\t5\tN\t5:0:0:0:0:0\t0:1:0:0:0:0\t1:0:0:0:0:0\n".
        "chr1\t6\tN\t6:0:0:0:0:0\t0:1:0:0:0:0\t1:0:0:0:0:0\n".
        "chr1\t7\tN\t7:0:0:0:0:0\t0:1:0:0:0:0\t1:0:0:0:0:0\n".
        "chr1\t8\tN\t8:0:0:0:0:0\t0:1:0:0:0:0\t1:0:0:0:0:0\n".
        "chr1\t9\tN\t9:0:0:0:0:0\t0:1:0:0:0:0\t10:0:0:0:0:0\n".
        "chr1\t10\tN\t10:0:0:0:0:0\t10:0:0:0:0:0\t10:0:0:0:0:0\n";
        $mc=get_max_coverage(\$file,"10%");
        is($mc->[0],9,"test maximum coverage method; maximum coverage is ok");
        is($mc->[1],1,"test maximum coverage script; maximum coverage is ok");
        is($mc->[2],1,"test maximum coverage script; maximum coverage is ok");
        
        $file=
        "chr1\t1\tN\t10:0:0:0:0:0\t0:1:0:0:0:0\t1:0:0:0:0:0\n".
        "chr1\t2\tN\t9:0:0:0:0:0\t0:1:0:0:0:0\t1:0:0:0:0:0\n".
        "chr1\t3\tN\t8:0:0:0:0:0\t0:1:0:0:0:0\t10:0:0:0:0:0\n".
        "chr1\t4\tN\t7:0:0:0:0:0\t0:10:0:0:0:0\t1:0:0:0:0:0\n".
        "chr1\t5\tN\t6:0:0:0:0:0\t0:1:0:0:0:0\t1:0:0:0:0:0\n".
        "chr1\t6\tN\t5:0:0:0:0:0\t0:10:0:0:0:0\t1:0:0:0:0:0\n".
        "chr1\t7\tN\t4:0:0:0:0:0\t0:1:0:0:0:0\t1:0:0:0:0:0\n".
        "chr1\t8\tN\t3:0:0:0:0:0\t0:1:0:0:0:0\t1:0:0:0:0:0\n".
        "chr1\t9\tN\t2:0:0:0:0:0\t0:1:0:0:0:0\t1:0:0:0:0:0\n".
        "chr1\t10\tN\t1:0:0:0:0:0\t1:0:0:0:0:0\t1:0:0:0:0:0\n";
        $mc=get_max_coverage(\$file,"10%");
        is($mc->[0],9,"test maximum coverage method; maximum coverage is ok");
        is($mc->[1],1,"test maximum coverage script; maximum coverage is ok");
        is($mc->[2],1,"test maximum coverage script; maximum coverage is ok");
        
        $file=
        "chr1\t1\tN\t10:0:0:0:0:0\t0:11:0:0:0:0\t21:0:0:0:0:0\n".
        "chr1\t2\tN\t9:0:0:0:0:0\t0:12:0:0:0:0\t22:0:0:0:0:0\n".
        "chr1\t3\tN\t8:0:0:0:0:0\t0:13:0:0:0:0\t23:0:0:0:0:0\n".
        "chr1\t4\tN\t7:0:0:0:0:0\t0:14:0:0:0:0\t24:0:0:0:0:0\n".
        "chr1\t5\tN\t6:0:0:0:0:0\t0:15:0:0:0:0\t25:0:0:0:0:0\n".
        "chr1\t6\tN\t5:0:0:0:0:0\t0:16:0:0:0:0\t26:0:0:0:0:0\n".
        "chr1\t7\tN\t4:0:0:0:0:0\t0:17:0:0:0:0\t27:0:0:0:0:0\n".
        "chr1\t8\tN\t3:0:0:0:0:0\t0:18:0:0:0:0\t28:0:0:0:0:0\n".
        "chr1\t9\tN\t2:0:0:0:0:0\t0:19:0:0:0:0\t29:0:0:0:0:0\n".
        "chr1\t10\tN\t1:0:0:0:0:0\t0:20:0:0:0:0\t30:0:0:0:0:0\n";
        $mc=get_max_coverage(\$file,"10%");
        is($mc->[0],9,"test maximum coverage method; maximum coverage is ok");
        is($mc->[1],19,"test maximum coverage script; maximum coverage is ok");
        is($mc->[2],29,"test maximum coverage script; maximum coverage is ok");
        
    }

    
}
1;