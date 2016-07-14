#!/usr/bin/perl
use strict;
use ProgramName;

my $MUMMIE=$ENV{"MUMMIE"};

my $name=ProgramName::get();
die "$name <in-dir> <model.hmm> <out-dir>\n" unless @ARGV==3;
my ($indir,$model,$outdir)=@ARGV;

my @files=`ls $indir`;
my $n=@files;
for(my $i=0 ; $i<$n ; ++$i) {
  my $file=$files[$i]; chomp $file;
  system("cp $indir/$file $outdir");
  open(OUT,">>$outdir/$file") || die;
  open(IN,"$MUMMIE/parse $model $indir/$file |") || die;
  <IN>; # skip header line
  while(<IN>) {
    print OUT $_;
  }
  close(IN);
  close(OUT);
}




