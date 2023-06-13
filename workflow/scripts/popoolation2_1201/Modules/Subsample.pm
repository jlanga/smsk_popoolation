{
package Subsample;
use strict;
use warnings;
use FindBin qw/$Bin/;
use lib "$Bin";
use Synchronized;

require Exporter;
our @ISA = qw(Exporter);
our @EXPORT     =qw(get_subsampler);




#withreplace, withoutreplace, fraction
sub get_subsampler
    {

        my $method=shift;
        my $targetcoverage=shift;
        $method=lc($method);
        
        my $subsampler;
        if($method eq "withreplace")
        {
            $subsampler=\&_withreplace;
        }
        elsif($method eq "withoutreplace")
        {
            $subsampler=\&_withoutreplace;
        }
        elsif($method eq "fraction")
        {
            $subsampler=\&_exactfraction;
        }
        else
        {
            die "unknown method for subsampling $method"
        }
        
        return sub
        {
            my $pileupentry=shift;
            return $subsampler->($pileupentry,$targetcoverage);
        }
    }



sub _withoutreplace
{
    my $p=shift;
    my $targetcov=shift;
    my($ac,$tc,$cc,$gc,$nc,$dc)=($p->{A},$p->{T},$p->{C},$p->{G},$p->{N},$p->{del});
    my $cov=$ac+$tc+$cc+$gc+$nc+$dc;
    return $p if $cov < $targetcov;
    
    

    
    
    my $temp="A"x$ac."T"x$tc."C"x$cc."G"x$gc."N"x$nc."*"x$dc;
    my @ar=split //,$temp;
    
    
    my($an,$tn,$cn,$gn,$nn,$dn)=(0,0,0,0,0,0);
    for my $i (1..$targetcov)
    {
        last unless @ar;
        my $idx=int(rand()*scalar(@ar));
        my $novel=splice(@ar,$idx,1);
        if($novel eq "A")
        {
            $an++;
        }
        elsif($novel eq "T")
        {
            $tn++;
        }
        elsif($novel eq "C")
        {
            $cn++;
        }
        elsif($novel eq "G")
        {
            $gn++;
        }
        elsif($novel eq "N")
        {
            $nn++;
        }
        elsif($novel eq "*")
        {
            $dn++;
        }
        else
        {
            die "invalid character for subsampled allele $novel";
        }
    }
    return{
        A=>$an,
        T=>$tn,
        C=>$cn,
        G=>$gn,
        N=>$nn,
        del=>$dn
    };
    
    
}


sub _withreplace
{
    my $p=shift;
    my $targetcov=shift;
    my($ac,$tc,$cc,$gc,$nc,$dc)=($p->{A},$p->{T},$p->{C},$p->{G},$p->{N},$p->{del});
    my $cov=$ac+$tc+$cc+$gc+$nc+$dc;
    return $p if $cov < $targetcov;
    # A T C G N del

    
    my($af,$tf,$cf,$gf,$nf,$df)=($ac/$cov , $tc/$cov , $cc/$cov , $gc/$cov , $nc/$cov , $dc/$cov);
    my $ab=$af;
    my $tb=$af+$tf;
    my $cb=$af+$tf+$cf;
    my $gb=$af+$tf+$cf+$gf;
    my $nb=$af+$tf+$cf+$gf+$nf;
    my $db=1;
    
    my($an,$tn,$cn,$gn,$nn,$dn)=(0,0,0,0,0,0);
    
    for my $i (1..$targetcov)
    {
        my $r=rand();
        if($r<$ab)
        {
            $an++;
        }
        elsif($r<$tb)
        {
            $tn++;
        }
        elsif($r<$cb)
        {
            $cn++;
        }
        elsif($r<$gb)
        {
            $gn++
        }
        elsif($r<$nb)
        {
            $nn++
        }
        elsif($r<$db)
        {
            $dn++
        }
        else
        {
            die "not valid random number $r must be 0<= random < 1"
        }
        
    }
    
    return
    {
        A=>$an,
        T=>$tn,
        C=>$cn,
        G=>$gn,
        N=>$nn,
        del=>$dn
    };
}
    
    sub _exactfraction
    {
        my $p=shift;
        my $targetcov=shift;
        my($ac,$tc,$cc,$gc,$nc,$dc)=($p->{A},$p->{T},$p->{C},$p->{G},$p->{N},$p->{del});
        my $cov=$ac+$tc+$cc+$gc+$nc+$dc;
        return $p if $cov < $targetcov;
        
        

        my($af,$tf,$cf,$gf,$nf,$df)=($ac/$cov , $tc/$cov , $cc/$cov , $gc/$cov , $nc/$cov , $dc/$cov);
        
        
        # ingenious, first iteratively increase the coverage until the target is reached, if coverage than is to high iteratively decrease it until targetcoveage is reached
        my $itercoverage=$targetcov;
        my($an,$tn,$cn,$gn,$nn,$dn)=(_intround($af*$itercoverage),_intround($tf*$itercoverage),_intround($cf*$itercoverage),
                                      _intround($gf*$itercoverage),_intround($nf*$itercoverage),_intround($df*$itercoverage));
        my $activecoverage=$an+$tn+$cn+$gn+$nn+$dn;
        while($activecoverage<$targetcov)
        {
            $itercoverage++;
            ($an,$tn,$cn,$gn,$nn,$dn)=(_intround($af*$itercoverage),_intround($tf*$itercoverage),_intround($cf*$itercoverage),
                                      _intround($gf*$itercoverage),_intround($nf*$itercoverage),_intround($df*$itercoverage));
            $activecoverage=$an+$tn+$cn+$gn+$nn+$dn;
        }
        while($activecoverage>$targetcov)
        {
            $itercoverage--;
            ($an,$tn,$cn,$gn,$nn,$dn)=(_intround($af*$itercoverage),_intround($tf*$itercoverage),_intround($cf*$itercoverage),
                                      _intround($gf*$itercoverage),_intround($nf*$itercoverage),_intround($df*$itercoverage));
            $activecoverage=$an+$tn+$cn+$gn+$nn+$dn;
        }
        
        my $tr=
         {
            A=>$an,
            T=>$tn,
            C=>$cn,
            G=>$gn,
            N=>$nn,
            del=>$dn
        };
        my $nuevocov=$tr->{A}+$tr->{T}+$tr->{C}+$tr->{G}+$tr->{N}+$tr->{del};
        die "calculated coverage is larger than the targetcoverage" if $nuevocov > $targetcov;
        
        # if coverage is still lower than the target fill it up with N
        $tr->{N}+=($targetcov-$nuevocov) if ($nuevocov<$targetcov);
        return $tr;
    }
    
    
    sub _intround
    {
        
        my $i=shift;            #2.5
        my $integer=int($i);    #2
        my $rest=$i-$integer;   #0.5
        if($rest<0.5)
        {
            return $integer;    #2
        }
        else{
            return $integer+1;  #3
        }
    }
}

1;