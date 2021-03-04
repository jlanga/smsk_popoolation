{
    package Test::TCMH;
    use strict;
    use warnings;
    use Statistics::R;
    use FindBin qw/$Bin/;
    use lib "$Bin/..";
    use Test;
    use MajorAlleles; # get the two major allele
    require Exporter;
    our @ISA = qw(Exporter);
    our @EXPORT  =qw(run_CMHTests);
    
    
    sub run_CMHTests
    {
        test_get2major_alleles();
    }
    


    
    
        
    sub test_get2major_alleles
    {
        
        #4	411770	T	0:25:0:6:0:0	0:28:0:2:0:0	0:24:0:5:0:0	0:36:0:7:0:0	0.4349568
        #4	411814	C	0:2:52:0:0:0	0:0:30:0:0:0	0:0:27:0:0:0	0:1:33:0:0:0	0.7454718
        #4	411917	A	38:1:0:0:0:0	40:0:0:0:0:0	16:0:0:0:0:0	11:2:0:0:0:0	0.9093106
        #4	412174	T	1:83:0:0:0:0	1:63:0:0:0:0	0:22:0:0:0:0	0:29:0:0:0:0	0.8465372
        #4	412187	G	2:0:0:84:0:0	3:0:0:53:0:0	0:0:0:25:0:0	2:0:0:30:0:0	0.2730484
        #4	412188	T	0:74:0:2:0:0	0:52:0:4:0:0	0:25:0:0:0:0	0:26:0:2:0:0	0.1683648
        #4	412199	A	76:0:0:0:0:0	54:0:0:0:0:0	29:0:0:0:0:0	33:0:0:2:0:0	0.5607681
        

        # 1
        my $data = [{samples=>[{A=>0,T=>25,C=>0,G=>6,N=>0,del=>0,eucov=>31},{A=>0,T=>28,C=>0,G=>2,N=>0,del=>0,eucov=>30},{A=>0,T=>24,C=>0,G=>5,N=>0,del=>0,eucov=>29},{A=>0,T=>36,C=>0,G=>7,N=>0,del=>0,eucov=>43}],ispuresnp=>1,refchar=>"T"}];
	my $samples = [{A=>0,T=>25,C=>0,G=>6,N=>0,del=>0,eucov=>31},{A=>0,T=>28,C=>0,G=>2,N=>0,del=>0,eucov=>30},{A=>0,T=>24,C=>0,G=>5,N=>0,del=>0,eucov=>29}];
	my ($major,$minor) = MajorAlleles::get_major_minor_alleles($samples);
        is($major,"T","Major allele one is correct");
        is($minor,"G","Major allele two is correct");
        
        # 2
        $data = [{samples=>[{A=>0,T=>2,C=>52,G=>0,N=>0,del=>0,eucov=>54},{A=>0,T=>0,C=>30,G=>0,N=>0,del=>0,eucov=>30},{A=>0,T=>0,C=>27,G=>0,N=>0,del=>0,eucov=>27},{A=>0,T=>1,C=>33,G=>0,N=>0,del=>0,eucov=>34}],ispuresnp=>1,refchar=>"C"}];
        $samples = [{A=>0,T=>2,C=>52,G=>0,N=>0,del=>0,eucov=>54},{A=>0,T=>0,C=>30,G=>0,N=>0,del=>0,eucov=>30},{A=>0,T=>0,C=>27,G=>0,N=>0,del=>0,eucov=>27},{A=>0,T=>1,C=>33,G=>0,N=>0,del=>0,eucov=>34}];
	($major,$minor) = MajorAlleles::get_major_minor_alleles($samples);
        is($major,"C","Major allele one is correct");
        is($minor,"T","Major allele two is correct");
        
        # 3
        $data = [{samples=>[{A=>38,T=>1,C=>0,G=>0,N=>0,del=>0,eucov=>39},{A=>40,T=>0,C=>0,G=>0,N=>0,del=>0,eucov=>40},{A=>16,T=>0,C=>0,G=>0,N=>0,del=>0,eucov=>16},{A=>11,T=>2,C=>0,G=>0,N=>0,del=>0,eucov=>13}],ispuresnp=>1,refchar=>"A"}];
        $samples = [{A=>38,T=>1,C=>0,G=>0,N=>0,del=>0,eucov=>39},{A=>40,T=>0,C=>0,G=>0,N=>0,del=>0,eucov=>40},{A=>16,T=>0,C=>0,G=>0,N=>0,del=>0,eucov=>16},{A=>11,T=>2,C=>0,G=>0,N=>0,del=>0,eucov=>13}];
	($major,$minor) = MajorAlleles::get_major_minor_alleles($samples);
        is($major,"A","Major allele one is correct");
        is($minor,"T","Major allele two is correct");
        
        # 4
        $data = [{samples=>[{A=>1,T=>83,C=>0,G=>0,N=>0,del=>0,eucov=>84},{A=>1,T=>63,C=>0,G=>0,N=>0,del=>0,eucov=>64},{A=>0,T=>22,C=>0,G=>0,N=>0,del=>0,eucov=>22},{A=>0,T=>29,C=>0,G=>0,N=>0,del=>0,eucov=>29}],ispuresnp=>1,refchar=>"T"}];
        $samples = [{A=>1,T=>83,C=>0,G=>0,N=>0,del=>0,eucov=>84},{A=>1,T=>63,C=>0,G=>0,N=>0,del=>0,eucov=>64},{A=>0,T=>22,C=>0,G=>0,N=>0,del=>0,eucov=>22},{A=>0,T=>29,C=>0,G=>0,N=>0,del=>0,eucov=>29}];
	($major,$minor) = MajorAlleles::get_major_minor_alleles($samples);
        is($major,"T","Major allele one is correct");
        is($minor,"A","Major allele two is correct");
        
        # 5
        $data = [{samples=>[{A=>2,T=>0,C=>0,G=>84,N=>0,del=>0,eucov=>86},{A=>3,T=>0,C=>0,G=>53,N=>0,del=>0,eucov=>56},{A=>0,T=>0,C=>0,G=>25,N=>0,del=>0,eucov=>25},{A=>2,T=>0,C=>0,G=>30,N=>0,del=>0,eucov=>32}],ispuresnp=>1,refchar=>"G"}];
	$samples = [{A=>2,T=>0,C=>0,G=>84,N=>0,del=>0,eucov=>86},{A=>3,T=>0,C=>0,G=>53,N=>0,del=>0,eucov=>56},{A=>0,T=>0,C=>0,G=>25,N=>0,del=>0,eucov=>25},{A=>2,T=>0,C=>0,G=>30,N=>0,del=>0,eucov=>32}];
        ($major,$minor) = MajorAlleles::get_major_minor_alleles($samples);
        is($major,"G","Major allele one is correct");
        is($minor,"A","Major allele two is correct");
        
        # 6
        $data = [{samples=>[{A=>0,T=>74,C=>0,G=>2,N=>0,del=>0,eucov=>76},{A=>0,T=>52,C=>0,G=>4,N=>0,del=>0,eucov=>56},{A=>0,T=>25,C=>0,G=>0,N=>0,del=>0,eucov=>25},{A=>0,T=>26,C=>0,G=>2,N=>0,del=>0,eucov=>28}],ispuresnp=>1,refchar=>"T"}];
        $samples = [{A=>0,T=>74,C=>0,G=>2,N=>0,del=>0,eucov=>76},{A=>0,T=>52,C=>0,G=>4,N=>0,del=>0,eucov=>56},{A=>0,T=>25,C=>0,G=>0,N=>0,del=>0,eucov=>25},{A=>0,T=>26,C=>0,G=>2,N=>0,del=>0,eucov=>28}];
	($major,$minor) = MajorAlleles::get_major_minor_alleles($samples);
        is($major,"T","Major allele one is correct");
        is($minor,"G","Major allele two is correct");
        
        # 7
        $data = [{samples=>[{A=>76,T=>0,C=>0,G=>0,N=>0,del=>0,eucov=>76},{A=>54,T=>0,C=>0,G=>0,N=>0,del=>0,eucov=>54},{A=>29,T=>0,C=>0,G=>0,N=>0,del=>0,eucov=>29},{A=>33,T=>0,C=>0,G=>2,N=>0,del=>0,eucov=>35}],ispuresnp=>1,refchar=>"A"}];
        $samples = [{A=>76,T=>0,C=>0,G=>0,N=>0,del=>0,eucov=>76},{A=>54,T=>0,C=>0,G=>0,N=>0,del=>0,eucov=>54},{A=>29,T=>0,C=>0,G=>0,N=>0,del=>0,eucov=>29},{A=>33,T=>0,C=>0,G=>2,N=>0,del=>0,eucov=>35}];
	($major,$minor) = MajorAlleles::get_major_minor_alleles($samples);
        is($major,"A","Major allele one is correct");
        is($minor,"G","Major allele two is correct");
        
        
        
    }
    
    
}
1;