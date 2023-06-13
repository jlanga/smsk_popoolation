{
    package SyncSlider;
    use strict;
    use warnings;
    use FindBin qw/$Bin/;
    use lib "$Bin";
    use List::Util qw[min max];
    use Synchronized;
    
    
    sub new
    {
        my $class=shift;
        my $file=shift;
        my $window=shift;
        my $step=shift;
        my $mincount=shift;
        my $mincov=shift;
        my $maxcov=shift;
        
        open my $fh,"<$file" or die "Could not open file handle";
        
        #get the parser of the synchronized file
        my $sp=get_sumsnp_synparser($mincount,$mincov,$maxcov);

        return bless {
            lower=>0,
            upper=>$window,
            window=>$window,
            step=>$step,
            file=>$file,
            fh=>$fh,
            sp=>$sp,
            curwin=>[],
            buffer=>[]
        },__PACKAGE__;
    }
    
    sub count_samples
    {
        my $self=shift;
        my $l=$self->_nextline();
        my $sp=$self->{sp};
        my $p=$sp->($l);
        my $c=scalar(@{$p->{samples}});
        $self->_bufferline($l);
        return $c;
    }
    
    sub nextWindow
    {
        my $self=shift;
        my $sp=$self->{sp}; # get the synchronized parser

        
        #get the current window, and the current chromosome
        my $curwin=$self->{curwin};
        
        my $curChr="";
        $curChr=$curwin->[0]{chr} if @$curwin;
        
        my $resetchr=0;
        
        # empty unnecessary entries
        EMPTY: while(@$curwin)
        {
            my $e=shift @$curwin;
            if($e->{pos}>$self->{lower})
            {
                unshift @$curwin, $e;
                last EMPTY;
            }
            
        }
        
        # fill with novel entries
        my $line;
        FILL:while($line=$self->_nextline)
        {
            my $e=$sp->($line);
            $curChr=$e->{chr} unless $curChr;
            
            
            if($e->{chr} eq $curChr && $e->{pos} <= $self->{upper})
            {
                push @$curwin,$e;
            }
            else
            {
                $resetchr=1 if $e->{chr} ne $curChr;
                $self->_bufferline($line);
                last FILL;
            }
        }
        
        return undef unless $curChr;
        
        
        my $toret=_annotateWindow($curwin,$curChr,$self->{lower},$self->{upper},$self->{window});
        
        if($resetchr or not defined($line))
        {
            # we transgressed the boundaries to the next chromosome
            # reset the windows and the current buffer
            $self->{lower}=0;
            $self->{upper}=$self->{window};
            $self->{curwin}=[];
        }
        else
        {
            # next time we will still be in the same chromosome
            # increase the upper and lower boundaries by the stepsize and set the current buffer
            $self->{upper}+=$self->{step};
            $self->{lower}+=$self->{step};
            $self->{curwin}=$curwin;
        }

        return $toret;
    }
    
    
    
    sub _nextline
    {
        my $self=shift;
        my $fh=$self->{fh};
        my $buffer=$self->{buffer};
        
        return shift @$buffer if @$buffer;
        return <$fh>;
    }
    
    sub _bufferline
    {
        my $self=shift;
        my $line=shift;
        push @{$self->{buffer}},$line;
    }
    
    
    
    sub _annotateWindow
    {
        my $curwin=shift;
        my $chr=shift;
        my $start=shift;
        my $end=shift;
        my $window=shift;
        
        my $snps=0;
        my $aboveCoverage=0;
        
        my $avmincov=0;
        foreach(@$curwin)
        {
            $snps++ if $_->{ispuresnp};
            if($_->{iscov})
            {
                $avmincov+=$_->{minsampcov};
                $aboveCoverage++;
            }
        }
        $avmincov/=$aboveCoverage if $aboveCoverage;

        return
        {
            chr=>$chr,
            start=>$start,
            end=>$end,
            middle=>int(($end+1+$start)/2),
            countpuresnp=>$snps,
            count_covered=>$aboveCoverage,
            window=>$window,
            avmincov=>$avmincov,
            data=>$curwin      
        };
    }
}
1;