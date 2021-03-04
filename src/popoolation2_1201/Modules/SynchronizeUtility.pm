{
    package SynchronizeUtility;
    use strict;
    use warnings;
    use FindBin qw/$Bin/;
    use lib "$Bin";
    use Synchronized;
    use Test;
    
    require Exporter;
    our @ISA = qw(Exporter);
    our @EXPORT =qw(format_parsed_pileup get_popcount_forsyncfile format_synchronized syncsample2string);
    
    
    sub syncsample2string
    {
        my $sync=shift;
        
        
        my $toret="A" x $sync->{A} . "T" x $sync->{T} . "C" x $sync->{C} . "G" x $sync->{G} . "-" x $sync->{del} . "N" x $sync->{N};
        return $toret; 
    }
    
    
    sub format_synchronized
    {
        my $sync=shift;
        
        my @temp=();
        # chr, pos, refchar, minsampcov, samples (A T C G N del)
        push @temp,$sync->{chr};
        push @temp,$sync->{pos};
        push @temp,$sync->{refchar};
        foreach my $s (@{$sync->{samples}})
        {
            
            push @temp,format_parsed_pileup($s);
        }
        
        my $toret=join("\t",@temp);
        return $toret;
    }
    
    sub format_parsed_pileup
    {
        my $pp=shift;
        return "$pp->{A}:$pp->{T}:$pp->{C}:$pp->{G}:$pp->{N}:$pp->{del}" if $pp;
        return "-";
    }
    
    sub get_popcount_forsyncfile
    {
        my $syncfile=shift;
        open my $ifh, "<", $syncfile or die "Could not open syncfile";
        my $pp=get_basic_syncparser();
        my $firstline=<$ifh>;
        chomp $firstline;
        my $pl=$pp->($firstline);
        my $count=scalar(@{$pl->{samples}});
        close $ifh;
        return $count;
        
        
    }
    


}

1;