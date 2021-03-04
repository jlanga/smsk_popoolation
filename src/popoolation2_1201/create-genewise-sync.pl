#!/usr/bin/env perl
{
use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use FindBin qw($RealBin);
use lib "$RealBin/Modules";
use Synchronized;

my $input="";
my $gtffile="";
my $output="";
my $help=0;
my $test=0;


GetOptions(
    "input=s"           =>\$input,
    "gtf=s"             =>\$gtffile,
    "output=s"          =>\$output,
    "test"              =>\$test,
    "help"              =>\$help
) or die "Invalid arguments";


pod2usage(-verbose=>2) if $help;
die "No test implemented for this scripts" if $test;
pod2usage(-msg=>"Could not find input file",-verbose=>1) unless -e $input;
pod2usage(-msg=>"Could not find gtf file",-verbose=>1) unless -e $gtffile;
pod2usage(-msg=>"Output file not provided",-verbose=>1) unless  $output;

my $paramfile=$output.".params";
open my $pfh, ">",$paramfile or die "Could not open $paramfile\n";
print $pfh "Using input\t$input\n";
print $pfh "Using output\t$output\n";
print $pfh "Using gtf\t$gtffile\n";
print $pfh "Using test\t$test\n";
print $pfh "Using help\t$help\n";
close $pfh;

my $syncparser=get_basic_syncparser();
my ($chrdec,$genehash)=Utility::read_gtf($gtffile);

print "Parsing sync file..\n";
open my $ifh, "<",$input or die "Could not open input file";
open my $ofh, ">",$output or die "Could not open output file";

while(my $line=<$ifh>)
{
    chomp $line;
    my($chr,$pos)=split /\t/,$line;
    my $genes=$chrdec->($chr,$pos);
    next unless $genes;
    
    my $sync=$syncparser->($line);
    
    foreach my $gene (@$genes)
    {
        my $gc=$genehash->{$gene};
        die "Error for gene $gene; position $sync->{pos} smaller than start $gc->{first}\n"if $sync->{pos} <   $gc->{first};
        die "Error for gene $gene; Position $sync->{pos} larger than end $gc->{last}\n"if $pos >   $gc->{last};
        push @{$gc->{sync}},$sync;
         
        # in case this is the last position of the gene, handle it immediately; print it         
        if($pos==$genehash->{$gene}{last})
        {
            Utility::handle_completed_gene($ofh,$genehash,$gene);
        }
    }
}

# handle the remaining genes/features
my @keys=keys(%$genehash);
for my $key (@keys)
{
    Utility::handle_completed_gene($ofh,$genehash,$key);
}


exit;
}

{
    use strict;
    use warnings;
    package Utility;
    use FindBin qw($RealBin);
    use lib "$RealBin/Modules";
    use SynchronizeUtility;
    
    sub handle_completed_gene
    {
        my $ofh=shift;
        my $genehash=shift;
        my $geneid=shift;
        
        my $temp=$genehash->{$geneid};
        #erase the entry from the hash -> very important memory issue
        delete($genehash->{$geneid});
        

        my $features=$temp->{features};
        my $first=$temp->{first};
        my $last=$temp->{last};
        my $strand=$temp->{strand};
        my $sync=$temp->{sync};
        return unless @$sync;

        my $binfeat={};
        foreach my $f (@$features)
        {
            my($start,$end)=($f->{start},$f->{end});
            ($start,$end)=($end,$start) if $start>$end;
            
            for my $i($start..$end)
            {
                $binfeat->{$i}=1;
            }
        }
        
        

        
        if($strand eq "+")
        {
            _write_plus_strand($ofh,$sync,$binfeat,$first,$geneid)
        }
        elsif($strand eq "-")
        {
            _write_minus_strand($ofh,$sync,$binfeat,$last,$geneid)
        }
        else
        {
            die "unknown strand $strand";
        }
    }
    
    sub _write_plus_strand
    {
        my $ofh=shift;
        my $syncs=shift;
        my $binfeat=shift;
        my $first=shift;
        my $geneid=shift;
        
        foreach my $s (@$syncs)
        {

            my $curpos=$s->{pos};
            my $nuevopos=0;
            for my $i ($first..$curpos)
            {
                $nuevopos++ if(exists($binfeat->{$i}));
            }
            
            my $nuevosyn=
            {
                chr=>$geneid,
                pos=>$nuevopos,
                refchar=>$s->{refchar},
                samples=>$s->{samples}
            };

            my $toprint=format_synchronized($nuevosyn);
            print $ofh  $toprint."\n";
        }
    }
    
    sub _write_minus_strand
    {
        my $ofh=shift;
        my $syncs=shift;
        my $binfeat=shift;
        my $last=shift;
        my $geneid=shift;
        $syncs=[reverse(@$syncs)];
        
        foreach my $s(@$syncs)
        {

            my $curpos=$s->{pos};
            my $nuevopos=0;
            for my $i ($curpos..$last)
            {
                $nuevopos++ if(exists($binfeat->{$i}));
            }

            my $nuevorefchar=$s->{refchar};
            $nuevorefchar=~tr/ATCGNatcgn/TAGCNtagcn/;
            
            my $nuevosamples=[];
            
            foreach my $samp (@{$s->{samples}})
            {
                my $ns={
                    A=>$samp->{T},
                    T=>$samp->{A},
                    C=>$samp->{G},
                    G=>$samp->{C},
                    N=>$samp->{N},
                    del=>$samp->{del}
                    };
                push @$nuevosamples, $ns;
            }
            
            my $nuevosyn=
            {
                chr=>$geneid,
                pos=>$nuevopos,
                refchar=>$nuevorefchar,
                samples=>$nuevosamples,
            };
            my $toprint=format_synchronized($nuevosyn);
            print $ofh $toprint."\n";
        }
    }
    
     sub _parsegtf
    {
        my $line=shift;
        my @a=split /\t/, $line;
        my $ref=$a[0];
        my $start=$a[3];
        my $end=$a[4];
        my $strand=$a[6];
        my $tfeat=$a[8];
            
        unless($ref or $start or $end or $tfeat)
        {
            die "the following line is not valid";
        }
        my $gene_id="";
        if($tfeat=~/gene_id "([^"]+)";/)
        {
            $gene_id=$1;
        }
        else
        {
            die "the following entry does not have a valid gene id: $line";
        }
        
        return
        {
            ref=>$ref,
            start=>$start,
            strand=>$strand,
            end=>$end,
            length=>$end-$start+1,
            gi=>$gene_id
        };
    }
    
    
    sub _getdefaultgenecoll
    {
        return
        {
            features=>[],
            length=>0,
            last=>0,
            strand=>undef,
            first=>1000000000000,
            sync=>[],
            covered=>0   
        }
        # last, length, covered, lines[], first, 
    }
    
    sub _getgeneint
    {
        my $geneid=shift;
        my $genemap=shift;
        my $lastcounter=shift;
        if(exists($genemap->{$geneid}))
        {
            return ($lastcounter,$genemap->{$geneid});
        }
        else
        {
            $lastcounter++;
            $genemap->{$geneid}=$lastcounter;
            return ($lastcounter,$genemap->{$geneid});
            
        }
    }
    
    sub _getDecodedGenemap
    {
        my $genemap=shift;
        
        my $decode=[];
        while(my($gene,$num)=each(%$genemap))
        {
            $decode->[$num]=$gene;
        }
        return $decode;
    }
    
    sub read_gtf
    {
        my $file=shift;
        open my $ifh,"<",$file or die "Could not open gtf-file";
        my $chrhash={};
        my $genecoll={};
        my $genemap={};
        my $lastcounter=0;

        
        print "Parsing gtf file..\n";
        while(my $line=<$ifh>)
        {
            chomp $line;
            my $ge=_parsegtf($line);
            
            my $gid=$ge->{gi};
            $genecoll->{$gid}=_getdefaultgenecoll unless exists($genecoll->{$gid});
            
            my $geneint;
            ($lastcounter,$geneint)=_getgeneint($gid,$genemap,$lastcounter);
            
            # update the chromosome hash
            my($start,$end,$ref)=($ge->{start},$ge->{end},$ge->{ref});
            push @{$genecoll->{$gid}{features}},
                {start=>$start,
                end=>$end};
            
            
            # set the new end if it is larger than the previous one
            $genecoll->{$gid}{last} = $end if $end > $genecoll->{$gid}{last};
            $genecoll->{$gid}{first} = $start if $start < $genecoll->{$gid}{first};
            $genecoll->{$gid}{strand} = $ge->{strand} unless($genecoll->{$gid}{strand});
            die "all features of a gene have to be on the same strand; Problem for $gid"  unless $genecoll->{$gid}{strand} eq $ge->{strand};
            
            #print "$ref $start $end\n";
            
            for(my $i=$start; $i<=$end; $i++)
            {
                if(exists($chrhash->{$ref}{$i}))
                {
                    my $ar=$chrhash->{$ref}{$i};
                    push @$ar,$geneint;
                    $ar=uniq($ar);
                    $chrhash->{$ref}{$i}=$ar;
                    
                }
                else
                {
                    $chrhash->{$ref}{$i}=[$geneint];
                }
            }
        }

        
        my $decodeGeneMap=_getDecodedGenemap($genemap);
################################################################################        
        # bless the beauty of a closure
        my $chrdecocer=sub {
            my $ref=shift;
            my $pos=shift;
            my $ta=$chrhash->{$ref}{$pos};
            return undef unless $ta;
            my $dec=[];
            for my $e (@$ta)
            {
                push @$dec,$decodeGeneMap->[$e];
            }
            return $dec;
        };
################################################################################
            
            
        #calculate the length of the features
        while(my($chr,$t)=each(%$chrhash))
        {
            while(my($pos,$genes)=each(%$t))
            {
                foreach my $g (@$genes)
                {
                    my $decg=$decodeGeneMap->[$g];
                    $genecoll->{$decg}{length}++;
                }
            }
        }
    
        return ($chrdecocer,$genecoll);
        #chr1 Twinscan  exon         501   650   .   +   .  gene_id "AB000381.000"; transcript_id "AB000381.000.1";
    }
    
    sub uniq
    {
        my $ar=shift;
        
        my $h={};
        foreach my $a (@$ar)
        {
            $h->{$a}=1;
        }
        return [keys(%$h)];
    }
    
}


=head1 NAME

perl create-genewise-sync.pl - Re-organizes a synchronized files by changing the coordinate system based on the provided features (e.g.: genes).

=head1 SYNOPSIS

perl create-genewise-sync.pl --input input.sync --gtf annotation.gtf --output genewise.sync

=head1 OPTIONS

=over 4

=item B<--input>

The input file in the synchronized format. Mandatory.

=item B<--gtf>

A gtf file as specified here: http://mblab.wustl.edu/GTF2.html

=item B<--output>

The output file in the synchronized format. Mandatory.

For more information use <--help>

=item B<--output>

The output file.  Mandatory.

=item B<--test>

Run the unit tests for this script. 

=item B<--help>

Display help for this script

=back

=head1 Details

=head2 Input synchronized


=head2 Input gtf

A gtf-file as described here http://mblab.wustl.edu/GTF2.html

 AB000381 Twinscan  exon         501   650   .   +   .  gene_id "AB000381.000"; transcript_id "AB000381.000.1";
 AB000381 Twinscan  CDS          501   650   .   +   2  gene_id "AB000381.000"; transcript_id "AB000381.000.1";
 AB000381 Twinscan  exon         700   800   .   +   .  gene_id "AB000381.000"; transcript_id "AB000381.000.1";
 AB000381 Twinscan  CDS          700   707   .   +   2  gene_id "AB000381.000"; transcript_id "AB000381.000.1";

The script is grouping all features having the same C<gene_id>! The feature field is not considered, so any feature may be used.
All features having the same 'gene_id' must also have the same strand. The tag C<transcript_id> will be ignored.
In case different features (eg. genes) are overlapping the respective entry is considered for every feature at a certain position once; 

=head2 Output



  
=cut
