#!/usr/bin/perl -w

use strict;
use warnings;
use Bio::Seq;
use Bio::SearchIO;
use List::MoreUtils qw(uniq);
use List::Util qw(min max);

# Parse HMMER output file
my $hmmer_file = "all_six_frames_smart.out";
my $in = Bio::SearchIO->new( -file => $hmmer_file, -format => 'hmmer' );

# Hashes to store start and end positions of hits
my %start_pos;
my %end_pos;

# Hash to store scaffold identifiers and their corresponding DNA sequence files
my %scaffold_files;

while (my $result = $in->next_result) {
    while (my $hit = $result->next_hit) {
        while (my $hsp = $hit->next_hsp) {
            if ($hsp->length('hit') > 65 && $hit->name() =~ /ig|I-set|Ig|V-set|IG/ && $hsp->evalue < 0.00001) {
                push(@{$start_pos{$result->query_name}}, $hsp->start('query'));
                push(@{$end_pos{$result->query_name}}, $hsp->end('query'));
            }
        }
    }
}

# Process hits and print results
print "strand\tscaffold\tscaf_start\tscaf_end\tchr_start\tchr_end\tIg-AA\tIg-nucl\n";
foreach my $key (sort keys %start_pos) {
    my $offset = min(@{$start_pos{$key}}) - 1;
    my $length = max(@{$end_pos{$key}}) - min(@{$start_pos{$key}});
    my $orf = $key;
    my $ig = substr($orf, $offset, $length);
    my $header = $key;
    $header =~ s/-/_/g;
    my ($scaf, $orf_start, $orf_end, $frame) = split("_", $header);
    my $adj_orf_start = $orf_start + $offset;
    my $adj_orf_end = $adj_orf_start + $length;

    # Infer scaffold file based on scaffold name
    my $scaffold_file = "$scaf.fasta";

    if (-e $scaffold_file) {
        open(my $fh, $scaffold_file) or die "Could not open $scaffold_file: $!";
        my $scaffold_seq = join('', <$fh>);
        close($fh);

        # Rest of your processing logic goes here...
        # Modify according to your specific requirements
        
        # Example: Print scaffold information
        print "Scaffold: $scaf\n";
        print "Scaffold Sequence: $scaffold_seq\n";
    } else {
        warn "Scaffold file $scaffold_file not found.";
    }
}

