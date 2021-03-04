{
package MajorAlleles;
use strict;
use warnings;
require Exporter;
our @ISA = qw(Exporter);
our @EXPORT     =qw(get_2major_alleles get_major_minor_alleles);

        sub get_major_minor_alleles
        {
                my $samples=shift;
                my ($ca,$ct,$cc,$cg,$cn,$cdel)=(0,0,0,0,0,0);
                foreach my $d (@$samples)
                {
                        $ca+= $d->{A};
                        $ct+= $d->{T};
                        $cc+= $d->{C};
                        $cg+= $d->{G};
                }
                
                my @als=({a=>"A",c=>$ca},{a=>"T",c=>$ct},{a=>"C",c=>$cc},{a=>"G",c=>$cg});
                @als=sort {$b->{c}<=>$a->{c}} @als;
                
                return ($als[0]{a},$als[1]{a});
        }



    sub get_2major_alleles {
        
                my $data=shift;

                my ($ca,$ct,$cc,$cg,$cn,$cdel)=(0,0,0,0,0,0);
                
                foreach my $d (@$data)
                {
                        next unless $d->{ispuresnp};
                        my $cosa=$d->{samples}; #countsamples
                        
                        for(my $i=0; $i<@$cosa; $i++) {
                                my $t=$cosa->[$i];
                                $ca+= $t->{A};
                                $ct+= $t->{T};
                                $cc+= $t->{C};
                                $cg+= $t->{G};
                                $cn+= $t->{N};
                                $cdel+= $t->{del};
                        }
                }
        
                my $act = $ca; #A count
                my $tct = $ct; #T count
                my $cct = $cc; #C count
                my $gct = $cg; #G count
        
        
                my @all_allele_count = ();
                my @vals = ();
                        
                $vals[0] = "A";
                $vals[1] = $act;
                push @all_allele_count, [ @vals ];
                @vals = ();
                        
                $vals[0] = "T";
                $vals[1] = $tct;
                push @all_allele_count, [ @vals ];
                @vals = ();
                        
                $vals[0] = "C";
                $vals[1] = $cct;
                push @all_allele_count, [ @vals ];
                @vals = ();
                        
                $vals[0] = "G";
                $vals[1] = $gct;
                push @all_allele_count, [ @vals ];
                @vals = ();
                        
                @all_allele_count = sort { $b->[1]<=> $a->[1] } @all_allele_count;
                    
           
                my $e = {
                        allele1=>$all_allele_count[0][0],
                        allele2=>$all_allele_count[1][0]
                };
                return $e;
        }
}
1;