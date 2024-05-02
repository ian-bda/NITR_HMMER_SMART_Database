#!/usr/bin/perl -w

use strict;
use warnings;
use Bio::Seq;
use Bio::SearchIO;
use List::MoreUtils qw(uniq);
use List::Util qw(min max);

# Configuration file containing scaffold identifiers and their corresponding DNA sequence files
my $config_file = "config.txt";

# Hash to store scaffold identifiers and their corresponding DNA sequence
my %scaffold_sequences;

# Read configuration file
open(CONFIG, $config_file) or die "Could not open $config_file: $!";
while (<CONFIG>) {
    chomp;
    my ($scaffold, $file) = split(/\s+/, $_, 2);
    open(my $fh, $file) or die "Could not open $file: $!";
    my $sequence = join('', <$fh>);
    close($fh);
    $scaffold_sequences{$scaffold} = $sequence;
}
close(CONFIG);

# Parse HMMER output file
my $hmmer_file = "all_six_frames_smart.out";
my $in = Bio::SearchIO->new( -file => $hmmer_file, -format => 'hmmer' );

# Hashes to store start and end positions of hits
my %start_pos;
my %end_pos;

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

    if (exists $scaffold_sequences{$scaf}) {
        my $scaffold_seq = $scaffold_sequences{$scaf};
        # Rest of your processing logic goes here...
        # Modify according to your specific requirements
    } else {
        warn "Scaffold $scaf not found in configuration file.";
    }
}
