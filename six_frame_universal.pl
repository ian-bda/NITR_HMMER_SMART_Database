#!/usr/bin/perl

### written by Ian Birchler de Allende on 03/29/23 based on code from Dustin Wcisel
### takes a scaffold fasta input and outputs all six frame translations with genomic coordinates embedded in the header
### subroutines from Beginning Perl for Bioinformatics, O'Reilly 2001

my $filename = $ARGV[0] or die "please provide a single sequence in fasta format\n";

# Extract filename without extension
my ($file_prefix) = $filename =~ /^(.+)\..+$/;

open(DNA, $filename);
while (<DNA>) { #this loop is necessary to remove any all blank lines, newlines, and the fasta header
	next if ($_ =~ /^>/); #ignores the line with the sequence name
	next if ($_ =~ /^\s*$/ ); #ignores any all blank lines
	chomp;
	push @dna, $_;
}

close (DNA);
my $dna = join ('', @dna);


###first frame
$orf = translate_frame($dna, 1);

$orig_peptide = $orf;
$orf =~ s/\*+/\*/g; #remove multiple stop codons and replace with a single

my @ORF = split(/\*/, $orf);

foreach (@ORF){
    if (length($_) > 70 && scalar(@count = $_ =~ /X/g)/length($_) < 0.1){ #removes short sequences and those with more than 10% X
        $orig_peptide =~ /$_/;
        print ">$file_prefix\_@--@+_Frame1\n$_\n";
    }
}

###second frame
$orf = translate_frame($dna, 2);

$orig_peptide = $orf;
$orf =~ s/\*+/\*/g; #remove multiple stop codons and replace with a single

my @ORF = split(/\*/, $orf);

foreach (@ORF){
    if (length($_) > 70 && scalar(@count = $_ =~ /X/g)/length($_) < 0.1){ #removes short sequences and those with more than 10% X
        $orig_peptide =~ /$_/;
        print ">$file_prefix\_@--@+_Frame2\n$_\n";
    }
}

###third frame
$orf = translate_frame($dna, 3);

$orig_peptide = $orf;
$orf =~ s/\*+/\*/g; #remove multiple stop codons and replace with a single

my @ORF = split(/\*/, $orf);

foreach (@ORF){
    if (length($_) > 70 && scalar(@count = $_ =~ /X/g)/length($_) < 0.1){ #removes short sequences and those with more than 10% X
        $orig_peptide =~ /$_/;
        print ">$file_prefix\_@--@+_Frame3\n$_\n";
    }
}

###fourth frame
$orf = translate_frame(revcom($dna), 1);

$orig_peptide = $orf;
$orf =~ s/\*+/\*/g; #remove multiple stop codons and replace with a single

my @ORF = split(/\*/, $orf);

foreach (@ORF){
    if (length($_) > 70 && scalar(@count = $_ =~ /X/g)/length($_) < 0.1){ #removes short sequences and those with more than 10% X
        $orig_peptide =~ /$_/;
        print ">$file_prefix\_@--@+_Frame4\n$_\n";
    }
}

###Fifth frame
$orf = translate_frame(revcom($dna), 2);

$orig_peptide = $orf;
$orf =~ s/\*+/\*/g; #remove multiple stop codons and replace with a single

my @ORF = split(/\*/, $orf);

foreach (@ORF){
    if (length($_) > 70 && scalar(@count = $_ =~ /X/g)/length($_) < 0.1){ #removes short sequences and those with more than 10% X
        $orig_peptide =~ /$_/;
        print ">$file_prefix\_@--@+_Frame5\n$_\n";
    }
}

###Sixth frame
$orf = translate_frame(revcom($dna), 3);

$orig_peptide = $orf;
$orf =~ s/\*+/\*/g; #remove multiple stop codons and replace with a single

my @ORF = split(/\*/, $orf);

foreach (@ORF){
    if (length($_) > 70 && scalar(@count = $_ =~ /X/g)/length($_) < 0.1){ #removes short sequences and those with more than 10% X
        $orig_peptide =~ /$_/;
        print ">$file_prefix\_@--@+_Frame6\n$_\n";
    }
}

sub codon2aa {
    my($codon) = @_;

    $codon = uc $codon;
 
    my(%genetic_code) = (
    
    'TCA' => 'S',    # Serine
    'TCC' => 'S',    # Serine
    'TCG' => 'S',    # Serine
    'TCT' => 'S',    # Serine
    'TTC' => 'F',    # Phenylalanine
    'TTT' => 'F',    # Phenylalanine
    'TTA' => 'L',    # Leucine
    'TTG' => 'L',    # Leucine
    'TAC' => 'Y',    # Tyrosine
    'TAT' => 'Y',    # Tyrosine
    'TAA' => '*',    # Stop
    'TAG' => '*',    # Stop
    'TGC' => 'C',    # Cysteine
    'TGT' => 'C',    # Cysteine
    'TGA' => '*',    # Stop
    'TGG' => 'W',    # Tryptophan
    'CTA' => 'L',    # Leucine
    'CTC' => 'L',    # Leucine
    'CTG' => 'L',    # Leucine
    'CTT' => 'L',    # Leucine
    'CCA' => 'P',    # Proline
    'CCC' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCT' => 'P',    # Proline
    'CAC' => 'H',    # Histidine
    'CAT' => 'H',    # Histidine
    'CAA' => 'Q',    # Glutamine
    'CAG' => 'Q',    # Glutamine
    'CGA' => 'R',    # Arginine
    'CGC' => 'R',    # Arginine
    'CGG' => 'R',    # Arginine
    'CGT' => 'R',    # Arginine
    'ATA' => 'I',    # Isoleucine
    'ATC' => 'I',    # Isoleucine
    'ATT' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine
    'ACA' => 'T',    # Threonine
    'ACC' => 'T',    # Threonine
    'ACG' => 'T',    # Threonine
    'ACT' => 'T',    # Threonine
    'AAC' => 'N',    # Asparagine
    'AAT' => 'N',    # Asparagine
    'AAA' => 'K',    # Lysine
    'AAG' => 'K',    # Lysine
    'AGC' => 'S',    # Serine
    'AGT' => 'S',    # Serine
    'AGA' => 'R',    # Arginine
    'AGG' => 'R',    # Arginine
    'GTA' => 'V',    # Valine
    'GTC' => 'V',    # Valine
    'GTG' => 'V',    # Valine
    'GTT' => 'V',    # Valine
    'GCA' => 'A',    # Alanine
    'GCC' => 'A',    # Alanine
    'GCG' => 'A',    # Alanine
    'GCT' => 'A',    # Alanine
    'GAC' => 'D',    # Aspartic Acid
    'GAT' => 'D',    # Aspartic Acid
    'GAA' => 'E',    # Glutamic Acid
    'GAG' => 'E',    # Glutamic Acid
    'GGA' => 'G',    # Glycine
    'GGC' => 'G',    # Glycine
    'GGG' => 'G',    # Glycine
    'GGT' => 'G',    # Glycine
    );

    if(exists $genetic_code{$codon}) {
        return $genetic_code{$codon};
    }else{
            return "X";
    }
}

sub translate_frame {

    my($seq, $start, $end) = @_;

    my $protein;

    # To make the subroutine easier to use, you won't need to specify
    #  the end point-it will just go to the end of the sequence
    #  by default.
    unless($end) {
        $end = length($seq);
    }

    # Finally, calculate and return the translation
        return dna2peptide ( substr ( $seq, $start - 1, $end -$start + 1) );
}

sub dna2peptide {

    my($dna) = @_;

    use strict;
    use warnings;

    # Initialize variables
    my $protein = '';

    # Translate each three-base codon to an amino acid, and append to a protein 
    for(my $i=0; $i < (length($dna) - 2) ; $i += 3) {
        $protein .= codon2aa( substr($dna,$i,3) );
    }

    return $protein;
}

sub revcom {

    my($dna) = @_;

    # First reverse the sequence
    my $revcom = reverse $dna;

    # Next, complement the sequence, dealing with upper and lower case
    # A->T, T->A, C->G, G->C
    $revcom =~ tr/ACGTacgt/TGCAtgca/;

    return $revcom;
}

