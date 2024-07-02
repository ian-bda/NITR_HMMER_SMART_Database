#!/usr/bin/perl -w

use strict;
use Bio::Seq;
use Bio::SearchIO;
use List::MoreUtils qw(uniq);
use List::Util qw(min max);
use Getopt::Long;

my %start_pos;
my %end_pos;

my $hmmer_file = "all_six_frames_smart.out";
my $all_six_frames = "all_six_frames.fa";
my @scaffold_files;

GetOptions(
    'hmmer_file=s'      => \$hmmer_file,
    'all_six_frames=s'  => \$all_six_frames,
    'scaffold=s'        => \@scaffold_files,
);

# Load HMMER output
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

# Load amino acid sequences
open (SIX_FRAMES, $all_six_frames);
my %headerhash;
my $header;

while (<SIX_FRAMES>){
    chomp;
    if ($_ =~ /^>/){
        $header = $_;
        $header =~ s/^>//g;
        next;
    }
    $headerhash{$header} = $_;
}
close (SIX_FRAMES);

# Load all DNA sequences into memory from input files
my %scaffold_sequences;
foreach my $fasta_file (@scaffold_files) {
    open(my $fh, '<', $fasta_file) or die "Cannot open $fasta_file: $!";
    my $scaffold_seq = '';
    while (<$fh>) {
        chomp;
        next if ($_ =~ /^>/); # ignores the line with the sequence name
        next if ($_ =~ /^\s*$/ ); # ignores any all blank lines
        $scaffold_seq .= $_;
    }
    close($fh);
    my ($scaffold_name) = $fasta_file =~ /^(.+)\.fasta$/;
    $scaffold_sequences{$scaffold_name} = $scaffold_seq;
}

my %helper = (
    'JANRMT010028174' => 0,
    'JANRMT010062807' => 0,
    'JANRMT010062854' => 0,
);

print "strand\tscaffold\tscaf_start\tscaf_end\tchr_start\tchr_end\tIg-AA\tIg-nucl\n";

foreach my $key (sort keys %start_pos) {
    my $offset = min(@{$start_pos{$key}}) - 1; # perl starts counting at 0
    my $length = max(@{$end_pos{$key}}) - min(@{$start_pos{$key}});
    my $orf = $headerhash{$key};
    my $ig = substr($orf, $offset, $length);
    my $header = $key;
    $header =~ s/-/_/g; # dash between numbers, change so split is easier
    my ($scaf, $orf_start, $orf_end, $frame) = split ("_", $header);
    my $adj_orf_start = $orf_start + $offset;
    my $adj_orf_end = $adj_orf_start + $length;

    if (exists $scaffold_sequences{$scaf}) {
        my $scaf_seq = $scaffold_sequences{$scaf};
        
        if ($frame =~ /Frame[1,2,3]/) { ## frames 1-3
            $frame =~ /(\d)/;
            my $dna_offset = ($adj_orf_start * 3) + $1 - 1;
            my $dna_length = ($adj_orf_end - $adj_orf_start) * 3;
            my $dna_end = $dna_offset + $dna_length;
            my $scaffold_seq = substr($scaf_seq, $dna_offset, $dna_length);
            my $chr_start = exists $helper{$scaf} ? $helper{$scaf} + $dna_offset : "NA";
            my $chr_end = exists $helper{$scaf} ? $helper{$scaf} + $dna_end : "NA";
            print "forward\t$scaf\t$dna_offset\t$dna_end\t$chr_start\t$chr_end\t$ig\t$scaffold_seq\n";
        }
        
        if ($frame =~ /Frame[4,5,6]/) { ## frames 4-6
            $frame =~ tr/456/123/;
            $frame =~ /(\d)/;
            my $revcom = reverse $scaf_seq;
            $revcom =~ tr/ACGTacgt/TGCAtgca/;
            my $dna_offset = ($adj_orf_start * 3) + $1 - 1;
            my $dna_length = ($adj_orf_end - $adj_orf_start) * 3;
            my $dna_end = $dna_offset + $dna_length;
            my $scaffold_seq = substr($revcom, $dna_offset, $dna_length);
            my $for_strand_start = length($scaf_seq) - $dna_offset;
            my $for_strand_end = $for_strand_start - $dna_length;
            my $chr_start = exists $helper{$scaf} ? $helper{$scaf} + $for_strand_start : "NA";
            my $chr_end = exists $helper{$scaf} ? $helper{$scaf} + $for_strand_end : "NA";
            print "reverse\t$scaf\t$for_strand_start\t$for_strand_end\t$chr_start\t$chr_end\t$ig\t$scaffold_seq\n";
        }
    }
}
