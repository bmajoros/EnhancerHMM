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
  print OUT "\%state\n";
  open(IN,"$MUMMIE/parse $model $indir/$file |") || die;
  <IN>; # skip header line
  while(<IN>) {
    chomp;
    my $state=$_+0;
    if($state==4) { $state=2 }
    elsif($state==5) { $state=3 }
    elsif($state==6) { $state=2 }
    elsif($state==7) { $state=1 }
    print OUT "$state\n";
  }
  close(IN);
  close(OUT);
}




