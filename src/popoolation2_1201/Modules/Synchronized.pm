{
    package Synchronized;
    use strict;
    use warnings;
    use FindBin qw/$Bin/;
    use lib "$Bin";
    use Test;
    use List::Util qw[min max];
    require Exporter;
    our @ISA = qw(Exporter);
    our @EXPORT     =qw(get_basic_syncparser get_sumsnp_synparser);
    
    
    #my $sp=get_basic_syncparser(),
    #my $line
    #my $pl=$sp->($line);
    # $pl->{chr}
    # $pl->{samples}[0]{A}
    
    
    sub get_basic_syncparser
    {
        return sub
        {
            my $line=shift;
            chomp $line;
            my @a=split /\s+/,$line;
            my $chr=shift @a;
            my $pos=shift @a;
            my $rc=shift @a;
            $rc=uc($rc);
            
            
            # Parsing
            my @samp;
            for(my $i=0; $i<@a; $i++)
            {
                my $col=$a[$i];
                my $e;
                
                if($col=~/-/)
                {
                    $e={index=>$i,eucov=>0,totcov=>0,A=>0,T=>0,C=>0,G=>0,N=>0,del=>0};
                }
                else
                {
                    my @parts=split /:/,$col;
                    die "failed parsing $col; does not have the correct number of entries" unless @parts ==6;
                    $e = {
                         A=>$parts[0],
                         T=>$parts[1],
                         C=>$parts[2],
                         G=>$parts[3],
                         N=>$parts[4],
                         del=>$parts[5]
                        };
                    
                    $e->{eucov}  = ($e->{A}+$e->{T}+$e->{C}+$e->{G});
                    $e->{totcov} = ($e->{eucov}+$e->{N}+$e->{del});
                    $e->{index}  = $i;
                }
                push @samp,$e;
            }
            
            my $minsampcov=min(map {$_->{eucov}} @samp);
            
            my $en={
            chr=>$chr,
            pos=>$pos,
            refchar=>$rc,
            minsampcov=>$minsampcov,
            samples=>\@samp
            };
            return $en;
        # chr, pos, refchar, minsampcov, samples (A T C G N del)
        }
    }
    
    sub get_sumsnp_synparser #mincount, mincoverage, maxcoverage
    {
        my $mincount=shift;
        my $mincoverage=shift;
        my $maxcoverage=shift;
        
        my $bp=get_basic_syncparser();
        
        return sub
        {
            my $line=shift;
            my $p=$bp->($line);
            
            
            # check if the position is a snp
            my @samp=@{$p->{samples}};
            my $issnp=0;
            my ($ca,$ct,$cc,$cg,$cn,$cdel)=(0,0,0,0,0,0);
            my $is_suficient_covered=1;
            
            
            for(my $i=0; $i<@samp; $i++)
            {
                my $s=$samp[$i];
                my $maxcov=$maxcoverage->[$i];
                my $iscov=1;
                $iscov=0 if($s->{eucov}<$mincoverage || $s->{eucov}>$maxcov);
                $s->{iscov}=$iscov;
                $ca+= $s->{A};
                $ct+= $s->{T};
                $cc+= $s->{C};
                $cg+= $s->{G};
                $cn+= $s->{N};
                $cdel+= $s->{del};
                
                # also reset the general appropriate coverage unless every entry has the correct coverage;
                $is_suficient_covered=0 unless $iscov;
                
            }
        
            my $allelecount=0;
            $allelecount++ if $ca >= $mincount;
            $allelecount++ if $ct >= $mincount;
            $allelecount++ if $cc >= $mincount;
            $allelecount++ if $cg >= $mincount;
            $issnp=1 if $allelecount > 1;          
        
        # unset the snp if not sufficiently covered
        $issnp=0 unless $is_suficient_covered; # no SNP will ever be allowed at a position which is not sufficiently covered in all data files!!
        
        # the SNP is tainted if there are any deletions at the given position
        my $taintedsnp = ($cdel >= $mincount)? 1:0;
        
        $p->{ignore}=$taintedsnp;
        $p->{iscov}=$is_suficient_covered;
        $p->{issnp}=$issnp;
        $p->{ispuresnp}=$p->{issnp};
        if($taintedsnp)
        {
            $p->{iscov}=0;
            $p->{ispuresnp}=0;
            $p->{issnp}=0;
        }
        
        # iscov, issnp, ispuresnps
        # chr, pos, refchar, minsampcov, samples (A T C G N del)
        
        return $p;
        };
    }
    
}
1;