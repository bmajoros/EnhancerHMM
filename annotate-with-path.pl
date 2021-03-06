#!/usr/bin/perl
use strict;
use ProgramName;

my $MUMMIE=$ENV{"MUMMIE"};

my $name=ProgramName::get();
die "$name <in-dir> <model.hmm> <0,1,2,3,2,1> <out-dir>\n" unless @ARGV==4;
my ($indir,$model,$classes,$outdir)=@ARGV;

my @STATE_TO_CLASS=split/,/,$classes;

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
    my $class=$STATE_TO_CLASS[$state];
    #if($state==4) { $state=2 }
    #elsif($state==5) { $state=1 }
#    elsif($state==6) { $state=2 }
#    elsif($state==7) { $state=1 }
    print OUT "$class\n";
  }
  close(IN);
  close(OUT);
}




