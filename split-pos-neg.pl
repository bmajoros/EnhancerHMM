#!/usr/bin/perl
use strict;

my $BASE="/home/bmajoros/GGR/p300";
my $POS_LIST="$BASE/positives.txt";
my $NEG_LIST="$BASE/negatives.txt";
my $POS_DIR="$BASE/pos";
my $NEG_DIR="$BASE/neg";
my $RAW="$BASE/raw";

copy($POS_LIST,$RAW,$POS_DIR);
copy($NEG_LIST,$RAW,$NEG_DIR);

#===============================================================
sub load {
  my ($infile,$indir,$outdir)=@_;
  open(IN,$infile) || die $infile;
  while(<IN>) {
    chomp; next unless(/\S/);
    my $file=$_;
    my $from="$indir/$file.t3.fastb";
    next unless -e $from;
    my $cmd="cp $from $outdir";
    System($cmd);
  }
  close(IN);
}
#===============================================================
sub System {
  my ($cmd)=@_;
  print "$cmd\n";
  # system($cmd);
}
#===============================================================
#===============================================================
#===============================================================
#===============================================================



