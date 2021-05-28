#!/usr/bin/perl
use strict;

my $BASE="/home/bmajoros/GGR/delta";
#my $JASPAR="/data/reddylab/software/meme_4.10.0/db/motif_databases/JASPAR/JASPAR_CORE_2016_vertebrates.meme";
my $JASPAR="/data/reddylab/projects/GGR/auxiliary/motifs/rsat_clustered/rsat_clustered_motifs.meme";
#my $OUTDIR="$BASE/mast-files";
my $OUTDIR="$BASE/meme-files";

my $header="MEME version 4

ALPHABET= ACGT

strands: + -

Background letter frequencies (from uniform background):
A 0.25000 C 0.25000 G 0.25000 T 0.25000 
";

open(IN,$JASPAR) || die "can't open file $JASPAR";
while(<IN>) {
  if(/^MOTIF\s+(\S+)\s+(\S+)/) {
    my ($id,$name)=($1,$2);
    my $outfile="$OUTDIR/$name.meme";
    open(OUT,">$outfile") || die "can't write to file $outfile";
    print OUT "$header\n";
    print OUT $_;
    while(<IN>) {
      last if(/URL/);
      print OUT
    }
    close(OUT);
  }
  elsif(/^MOTIF\s+(\S+)/) {
    my ($id)=($1);
    my $name=$id;
    my $outfile="$OUTDIR/$name.meme";
    open(OUT,">$outfile") || die "can't write to file $outfile";
    print OUT "$header\n";
    print OUT $_;
    $_=<IN>; print OUT $_;
    while(<IN>) {
      print OUT;
      last unless(/\S/);
    }
    close(OUT);
  }
}
close(IN);




