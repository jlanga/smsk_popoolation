{
    package Test::TFstCalculation;
    use strict;
    use warnings;
    use FindBin qw/$Bin/;
    use lib "$Bin/..";
    use Test;
    use FstAsymptUnbiased qw/get_asymptunbiased_fstcalculator reorganizeSNPs calculate_nk_dk/;
    require Exporter;
    our @ISA = qw(Exporter);
    our @EXPORT  =qw(run_FstCalculationTests);
    
    
    sub run_FstCalculationTests
    {
        test_reorganize_snp();
        test_calculate_nk_dk();
    }
    
    
    sub test_reorganize_snp
    {
        my $re;
        my $snpar;
        my $t;
        
        $re=reorganizeSNPs([{samples=>[{A=>2,T=>2,C=>0,G=>0,eucov=>4},{A=>3,T=>3,C=>0,G=>0,eucov=>6}]}],2);
        $t=$re->[0]{"0:1"};
        
        is($t->{a1},2,"Karlsson Fst calculation: SNP reorganization correct");
        is($t->{a2},3,"Karlsson Fst calculation: SNP reorganization correct");
        is($t->{n1},4,"Karlsson Fst calculation: SNP reorganization correct");
        is($t->{n2},6,"Karlsson Fst calculation: SNP reorganization correct");
        
        $re=reorganizeSNPs([{samples=>[{A=>2,T=>2,C=>2,G=>2,eucov=>8},{A=>3,T=>3,C=>0,G=>0,eucov=>6}]}],2);
        $t=$re->[0]{"0:1"};
        
        is($t->{a1},2,"Karlsson Fst calculation: SNP reorganization correct");
        is($t->{a2},3,"Karlsson Fst calculation: SNP reorganization correct");
        is($t->{n1},4,"Karlsson Fst calculation: SNP reorganization correct");
        is($t->{n2},6,"Karlsson Fst calculation: SNP reorganization correct");
        
        
        $re=reorganizeSNPs([{samples=>[{A=>0,T=>3,C=>1,G=>0,eucov=>4},{A=>3,T=>3,C=>0,G=>0,eucov=>6}]}],2);
        $t=$re->[0]{"0:1"};
        
        is($t->{a1},3,"Karlsson Fst calculation: SNP reorganization correct");
        is($t->{a2},3,"Karlsson Fst calculation: SNP reorganization correct");
        is($t->{n1},3,"Karlsson Fst calculation: SNP reorganization correct");
        is($t->{n2},6,"Karlsson Fst calculation: SNP reorganization correct");
        
        $re=reorganizeSNPs([{samples=>[{A=>2,T=>2,C=>1,G=>1,eucov=>6},{A=>3,T=>3,C=>0,G=>0,eucov=>6}, {A=>0,T=>0,C=>4,G=>4,eucov=>8} ]},],3);
        $t=$re->[0]{"0:1"};
        is($t->{a1},2,"Karlsson Fst calculation: SNP reorganization correct");
        is($t->{a2},3,"Karlsson Fst calculation: SNP reorganization correct");
        is($t->{n1},4,"Karlsson Fst calculation: SNP reorganization correct");
        is($t->{n2},6,"Karlsson Fst calculation: SNP reorganization correct");
        
        $t=$re->[0]{"0:2"};
        is($t->{a1},1,"Karlsson Fst calculation: SNP reorganization correct");
        is($t->{a2},4,"Karlsson Fst calculation: SNP reorganization correct");
        is($t->{n1},2,"Karlsson Fst calculation: SNP reorganization correct");
        is($t->{n2},8,"Karlsson Fst calculation: SNP reorganization correct");
        
        $t=$re->[0]{"1:2"};
        is($t->{a1},3,"Karlsson Fst calculation: SNP reorganization correct");
        is($t->{a2},0,"Karlsson Fst calculation: SNP reorganization correct");
        is($t->{n1},6,"Karlsson Fst calculation: SNP reorganization correct");
        is($t->{n2},0,"Karlsson Fst calculation: SNP reorganization correct");
    }
    
    sub test_calculate_nk_dk
    {
        my ($nk,$dk);
        
        ($nk,$dk)=calculate_nk_dk({a1=>2,n1=>4,a2=>2,n2=>4});
        ok(abs($nk-(-0.16666))<0.0001,"Karlsson Fst calculation: Nk calculated correctly");
        ok(abs($dk-(0.5))<0.0001,"Karlsson Fst calculation: Dk calculated correctly");
        
        ($nk,$dk)=calculate_nk_dk({a1=>4,n1=>8,a2=>4,n2=>8});
        ok(abs($nk-(-0.0714286))<0.0001,"Karlsson Fst calculation: Nk calculated correctly");
        ok(abs($dk-(0.5))<0.0001,"KarlssonFst calculation: Dk calculated correctly");

        ($nk,$dk)=calculate_nk_dk({a1=>200,n1=>400,a2=>4,n2=>8});
        ok(abs($nk-(-0.0363409))<0.0001,"Karlsson Fst calculation: Nk calculated correctly");
        ok(abs($dk-(0.5))<0.0001,"Karlsson Fst calculation: Dk calculated correctly");

        ($nk,$dk)=calculate_nk_dk({a1=>5,n1=>10,a2=>1,n2=>10});
        ok(abs($nk-(0.12222222))<0.0001,"Karlsson Fst calculation: Nk calculated correctly");
        ok(abs($dk-(0.5))<0.0001,"Karlsson Fst calculation: Dk calculated correctly");

        ($nk,$dk)=calculate_nk_dk({a1=>1,n1=>10,a2=>2,n2=>100});
        ok(abs($nk-(-0.003798))<0.0001,"Karlsson Fst calculation: Nk calculated correctly");
        ok(abs($dk-(0.116))<0.0001,"Karlsson Fst calculation: Dk calculated correctly");
    }
    
    
}
1;