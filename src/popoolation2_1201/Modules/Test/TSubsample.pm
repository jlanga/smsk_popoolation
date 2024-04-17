{
    package Test::TSubsample;
    use strict;
    use warnings;
    use FindBin qw/$Bin/;
    use lib "$Bin/..";
    use Test;
    use Subsample;
    require Exporter;
    our @ISA = qw(Exporter);
    our @EXPORT  =qw(run_SubsampleTests);
    
    
    sub run_SubsampleTests
    {
        test_withreplacement();
        test_withoutreplacement();
        test_fraction();
    }
    
    
    sub test_withreplacement
    {
       my $ss;
       my $r;
       
       $ss=get_subsampler("withreplace",10);
       $r=$ss->({A=>20,T=>0,C=>0,G=>0,N=>0,del=>0});
       is($r->{A},10,"Test subsampler with replacment, A count is correct");
       is($r->{T},0, "Test subsampler with replacment, T count is correct");
       is($r->{C},0, "Test subsampler with replacment, C count is correct");
       is($r->{G},0, "Test subsampler with replacment, G count is correct");
        $r=$ss->({A=>0,T=>20,C=>0,G=>0,N=>0,del=>0});
       is($r->{A},0,"Test subsampler with replacment, A count is correct");
       is($r->{T},10, "Test subsampler with replacment, T count is correct");
       is($r->{C},0, "Test subsampler with replacment, C count is correct");
       is($r->{G},0, "Test subsampler with replacment, G count is correct");
        $r=$ss->({A=>0,T=>0,C=>20,G=>0,N=>0,del=>0});
       is($r->{A},0,"Test subsampler with replacment, A count is correct");
       is($r->{T},0, "Test subsampler with replacment, T count is correct");
       is($r->{C},10, "Test subsampler with replacment, C count is correct");
       is($r->{G},0, "Test subsampler with replacment, G count is correct");
        $r=$ss->({A=>0,T=>0,C=>0,G=>20,N=>0,del=>0});
       is($r->{A},0,"Test subsampler with replacment, A count is correct");
       is($r->{T},0, "Test subsampler with replacment, T count is correct");
       is($r->{C},0, "Test subsampler with replacment, C count is correct");
       is($r->{G},10, "Test subsampler with replacment, G count is correct");
       
       # coverage to small
        $r=$ss->({A=>2,T=>2,C=>2,G=>2,N=>0,del=>0});
       is($r->{A},2,"Test subsampler with replacment, A count is correct");
       is($r->{T},2, "Test subsampler with replacment, T count is correct");
       is($r->{C},2, "Test subsampler with replacment, C count is correct");
       is($r->{G},2, "Test subsampler with replacment, G count is correct");
       
       
       
    }
    
    sub test_withoutreplacement
    {
        my $ss;
       my $r;
        $ss=get_subsampler("withoutreplace",10);
       $r=$ss->({A=>20,T=>0,C=>0,G=>0,N=>0,del=>0});
       is($r->{A},10,"Test subsampler without replacment, A count is correct");
       is($r->{T},0, "Test subsampler without replacment, T count is correct");
       is($r->{C},0, "Test subsampler without replacment, C count is correct");
       is($r->{G},0, "Test subsampler without replacment, G count is correct");
        $r=$ss->({A=>0,T=>20,C=>0,G=>0,N=>0,del=>0});
       is($r->{A},0,"Test subsampler without replacment, A count is correct");
       is($r->{T},10, "Test subsampler without replacment, T count is correct");
       is($r->{C},0, "Test subsampler without replacment, C count is correct");
       is($r->{G},0, "Test subsampler without replacment, G count is correct");
        $r=$ss->({A=>0,T=>0,C=>20,G=>0,N=>0,del=>0});
       is($r->{A},0,"Test subsampler without replacment, A count is correct");
       is($r->{T},0, "Test subsampler without replacment, T count is correct");
       is($r->{C},10, "Test subsampler without replacment, C count is correct");
       is($r->{G},0, "Test subsampler without replacment, G count is correct");
        $r=$ss->({A=>0,T=>0,C=>0,G=>20,N=>0,del=>0});
       is($r->{A},0,"Test subsampler without replacment, A count is correct");
       is($r->{T},0, "Test subsampler without replacment, T count is correct");
       is($r->{C},0, "Test subsampler without replacment, C count is correct");
       is($r->{G},10, "Test subsampler without replacment, G count is correct");
       
       # coverage to small
        $r=$ss->({A=>2,T=>2,C=>2,G=>2,N=>0,del=>0});
       is($r->{A},2,"Test subsampler without replacment, A count is correct");
       is($r->{T},2, "Test subsampler without replacment, T count is correct");
       is($r->{C},2, "Test subsampler without replacment, C count is correct");
       is($r->{G},2, "Test subsampler without replacment, G count is correct");
       
    }
    
    sub test_fraction
    {
     
        my $ss;
        my $r;
        $ss=get_subsampler("fraction",10);
       $r=$ss->({A=>20,T=>0,C=>0,G=>0,N=>0,del=>0});
       is($r->{A},10,"Test subsampler fraction, A count is correct");
       is($r->{T},0, "Test subsampler fraction, T count is correct");
       is($r->{C},0, "Test subsampler fraction, C count is correct");
       is($r->{G},0, "Test subsampler fraction, G count is correct");
        $r=$ss->({A=>0,T=>20,C=>0,G=>0,N=>0,del=>0});
       is($r->{A},0,"Test subsampler fraction, A count is correct");
       is($r->{T},10, "Test subsampler fraction, T count is correct");
       is($r->{C},0, "Test subsampler fraction, C count is correct");
       is($r->{G},0, "Test subsampler fraction, G count is correct");
        $r=$ss->({A=>0,T=>0,C=>20,G=>0,N=>0,del=>0});
       is($r->{A},0,"Test subsampler fraction, A count is correct");
       is($r->{T},0, "Test subsampler fraction, T count is correct");
       is($r->{C},10, "Test subsampler fraction, C count is correct");
       is($r->{G},0, "Test subsampler fraction, G count is correct");
        $r=$ss->({A=>0,T=>0,C=>0,G=>20,N=>0,del=>0});
       is($r->{A},0,"Test subsampler fraction, A count is correct");
       is($r->{T},0, "Test subsampler fraction, T count is correct");
       is($r->{C},0, "Test subsampler fraction, C count is correct");
       is($r->{G},10, "Test subsampler fraction, G count is correct");
       
       # coverage to small
        $r=$ss->({A=>2,T=>2,C=>2,G=>2,N=>0,del=>0});
       is($r->{A},2,"Test subsampler fraction, A count is correct");
       is($r->{T},2, "Test subsampler fraction, T count is correct");
       is($r->{C},2, "Test subsampler fraction, C count is correct");
       is($r->{G},2, "Test subsampler fraction, G count is correct");
    }
    
    
}
1;