#!/usr/bin/perl
use strict;

my $MAX_FILES=1000;
my $BASE="/home/bmajoros/GGR/p300";
my $POS_LIST="$BASE/positives.txt";
my $NEG_LIST="$BASE/negatives.txt";
my $POS_DIR="$BASE/pos";
my $NEG_DIR="$BASE/neg";
my $RAW="$BASE/raw";

copy($POS_LIST,$RAW,$POS_DIR);
copy($NEG_LIST,$RAW,$NEG_DIR);

#===============================================================
sub copy {
  my ($infile,$indir,$outdir)=@_;
  open(IN,$infile) || die $infile;
  my $copied=0;
  while(<IN>) {
    chomp; next unless(/\S/);
    my $file=$_;
    my $from="$indir/$file.t3.fastb";
    next unless -e $from;
    my $cmd="cp $from $outdir";
    System($cmd);
    ++$copied;
    last if $copied>=$MAX_FILES;
  }
  close(IN);
}
#===============================================================
sub System {
  my ($cmd)=@_;
  print "$cmd\n";
  system($cmd);
}
#===============================================================
#===============================================================
#===============================================================
#===============================================================



