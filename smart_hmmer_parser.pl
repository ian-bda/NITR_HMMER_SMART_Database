#!/usr/bin/perl -w

use Bio::Seq;
use Bio::SearchIO;
use List::MoreUtils qw(uniq);
use strict;
use List::Util qw(min max);

my %start_pos;
my %end_pos;

my $hmmer_file = "all_six_frames_smart.out";
my $all_six_frames = "all_six_frames.fa";

my $in = Bio::SearchIO->new( -file => $hmmer_file, -format => 'hmmer' );
while ( my $result = $in -> next_result ) {
    while( my $hit = $result->next_hit ) {
        while( my $hsp = $hit->next_hsp ) { 
		    if ($hsp->length('hit') > 65 && $hit->name() =~ /ig|I-set|Ig|V-set|IG/ &&  $hsp->evalue < 0.00001){
	            push ( @{$start_pos{$result->query_name}}, $hsp->start('query') );
	            push ( @{$end_pos{$result->query_name}}, $hsp->end('query') );
            }
        }
    }
}

open (SIX_FRAMES, $all_six_frames);

my %headerhash;
my $header;

while (<SIX_FRAMES>){ #need this to get the full length amino acid sequences to substring later
    chomp;
    if ($_ =~ /^>/){
        $header = $_;
        $header =~ s/^>//g;
        next;
    }
    $headerhash{$header} = $_;
}

close (SIX_FRAMES);

###load all the DNA sequence into memory
my $JANRMT010028174_file = "JANRMT010028174.fasta";
my $JANRMT010062807_file = "JANRMT010062807.fasta";
my $JANRMT010062854_file = "JANRMT010062854.fasta";

my @dna1;
open(DNA1, $JANRMT010028174_file);
while (<DNA1>) { #necessary to remove any all blank lines, newlines, and the fasta header
	next if ($_ =~ /^>/); #ignores the line with the sequence name
	next if ($_ =~ /^\s*$/ ); #ignores any all blank lines
	chomp;
	push @dna1, $_;
}

close (DNA1);
my $scaf_JANRMT010028174 = join ('', @dna1);

my @dna2;
open(DNA2, $JANRMT010062807_file);
while (<DNA2>) { #necessary to remove any all blank lines, newlines, and the fasta header
	next if ($_ =~ /^>/); #ignores the line with the sequence name
	next if ($_ =~ /^\s*$/ ); #ignores any all blank lines
	chomp;
	push @dna2, $_;
}

close (DNA2);
my $scaf_JANRMT010062807 = join ('', @dna2);

my @dna3;
open(DNA3, $JANRMT010062854_file);
while (<DNA3>) { #necessary to remove any all blank lines, newlines, and the fasta header
	next if ($_ =~ /^>/); #ignores the line with the sequence name
	next if ($_ =~ /^\s*$/ ); #ignores any all blank lines
	chomp;
	push @dna3, $_;
}

close (DNA3);
my $scaf_JANRMT010062854 = join ('', @dna3);
###finished getting DNA sequences for all three scaffolds

my $helper101 = 0;
my $helper100 = 0;
my $helper111 = 0;

print "strand\tscaffold\tscaf_start\tscaf_end\tchr_start\tchr_end\tIg-AA\tIg-nucl\n";

foreach my $key (sort keys %start_pos){
    my $offset = min(@{$start_pos{$key}})-1; #perl starts counting at 0
    my $length = max(@{$end_pos{$key}}) - min(@{$start_pos{$key}});
    my $orf = $headerhash{$key};
    my $ig = substr($orf, $offset, $length);
    my $header = $key;
    $header =~ s/-/_/g; #dash between numbers, change so split is easier
    my ($scaf, $orf_start, $orf_end, $frame) = split ("_", $header);
    my $adj_orf_start = $orf_start+$offset;
    my $adj_orf_end = $adj_orf_start+$length;

    if ($scaf =~ /^JANRMT010028174/ && $frame =~ /Frame[1,2,3]/){ ## frames 1-3
        $frame =~ /(\d)/;
        my $dna_offset = ($adj_orf_start*3)+$1-1;
        my $dna_length = ($adj_orf_end-$adj_orf_start)*3;
        my $dna_end = $dna_offset+$dna_length;
        my $scaffold_seq = substr($scaf_JANRMT010028174, $dna_offset, $dna_length);
        my $chr_start = $helper111 + $dna_offset;
        my $chr_end = $helper111 + $dna_end;
        print "forward\t$scaf\t$dna_offset\t$dna_end\t$chr_start\t$chr_end\t$ig\t$scaffold_seq\n";
    }
    if ($scaf =~ /^JANRMT010028174/ && $frame =~ /Frame[4,5,6]/){ ## frames 4-6
        $frame =~ tr/456/123/;
        $frame =~ /(\d)/;
        my $revcom = reverse $scaf_JANRMT010028174;
        $revcom =~ tr/ACGTacgt/TGCAtgca/;
        my $dna_offset = ($adj_orf_start*3)+$1-1;
        my $dna_length = ($adj_orf_end-$adj_orf_start)*3;
        my $dna_end = $dna_offset+$dna_length;
        my $scaffold_seq = substr($revcom, $dna_offset, $dna_length);
        my $for_strand_start = length($scaf_JANRMT010028174)-$dna_offset;
        my $for_strand_end = $for_strand_start-$dna_length;
        my $chr_start = $helper111 + $for_strand_start;
        my $chr_end = $helper111 + $for_strand_end;
        print "reverse\t$scaf\t$for_strand_start\t$for_strand_end\t$chr_start\t$chr_end\t$ig\t$scaffold_seq\n";
    }
    if ($scaf =~ /^JANRMT010062807/ && $frame =~ /Frame[1,2,3]/){ ## frames 1-3
        $frame =~ /(\d)/;
        my $dna_offset = ($adj_orf_start*3)+$1-1;
        my $dna_length = ($adj_orf_end-$adj_orf_start)*3;
        my $dna_end = $dna_offset+$dna_length;
        my $scaffold_seq = substr($scaf_JANRMT010062807, $dna_offset, $dna_length);
        my $chr_start = $helper101 + $dna_offset;
        my $chr_end = $helper101 + $dna_end;
        print "forward\t$scaf\t$dna_offset\t$dna_end\t$chr_start\t$chr_end\t$ig\t$scaffold_seq\n";
    }
    if ($scaf =~ /^JANRMT010062807/ && $frame =~ /Frame[4,5,6]/){ ## frames 4-6
        $frame =~ tr/456/123/;
        $frame =~ /(\d)/;
        my $revcom = reverse $scaf_JANRMT010062807;
        $revcom =~ tr/ACGTacgt/TGCAtgca/;
        my $dna_offset = ($adj_orf_start*3)+$1-1;
        my $dna_length = ($adj_orf_end-$adj_orf_start)*3;
        my $dna_end = $dna_offset+$dna_length;
        my $scaffold_seq = substr($revcom, $dna_offset, $dna_length);
        my $for_strand_start = length($scaf_JANRMT010062807)-$dna_offset;
        my $for_strand_end = $for_strand_start-$dna_length;
        my $chr_start = $helper101 + $for_strand_start;
        my $chr_end = $helper101 + $for_strand_end;
        print "reverse\t$scaf\t$for_strand_start\t$for_strand_end\t$chr_start\t$chr_end\t$ig\t$scaffold_seq\n";
    }
    if ($scaf =~ /^JANRMT010062854/ && $frame =~ /Frame[1,2,3]/){ ## frames 1-3
        $frame =~ /(\d)/;
        my $dna_offset = ($adj_orf_start*3)+$1-1;
        my $dna_length = ($adj_orf_end-$adj_orf_start)*3;
        my $dna_end = $dna_offset+$dna_length;
        my $scaffold_seq = substr($scaf_JANRMT010062854, $dna_offset, $dna_length);
        my $chr_start = $helper100 + $dna_offset;
        my $chr_end = $helper100 + $dna_end;
        print "forward\t$scaf\t$dna_offset\t$dna_end\t$chr_start\t$chr_end\t$ig\t$scaffold_seq\n";
    }
    if ($scaf =~ /^JANRMT010062854/ && $frame =~ /Frame[4,5,6]/){ ## frames 4-6
        $frame =~ tr/456/123/;
        $frame =~ /(\d)/;
        my $revcom = reverse $scaf_JANRMT010062854;
        $revcom =~ tr/ACGTacgt/TGCAtgca/;
        my $dna_offset = ($adj_orf_start*3)+$1-1;
        my $dna_length = ($adj_orf_end-$adj_orf_start)*3;
        my $dna_end = $dna_offset+$dna_length;
        my $scaffold_seq = substr($revcom, $dna_offset, $dna_length);
        my $for_strand_start = length($scaf_JANRMT010062854)-$dna_offset;
        my $for_strand_end = $for_strand_start-$dna_length;
        my $chr_start = $helper100 + $for_strand_start;
        my $chr_end = $helper100 + $for_strand_end;
        print "reverse\t$scaf\t$for_strand_start\t$for_strand_end\t$chr_start\t$chr_end\t$ig\t$scaffold_seq\n";
    }
}

