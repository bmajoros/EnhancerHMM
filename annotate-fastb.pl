#!/usr/bin/perl
use strict;
use ProgramName;

my $MUMMIE=$ENV{"MUMMIE"};

my $name=ProgramName::get();
die "$name <in.fastb> <model.hmm> <out.fastb>\n" unless @ARGV==3;
my ($infile,$model,$outfile)=@ARGV;

system("cp $infile $outfile");
open(OUT,">>$outfile") || die;
print OUT "\%state\n";
open(IN,"$MUMMIE/parse $model $infile |") || die;
<IN>; # skip header line
while(<IN>) {
  chomp;
  my $state=$_+0;
  if($state==4) { $state=2 }
  elsif($state==5) { $state=1 }
  #elsif($state==6) { $state=2 }
  #elsif($state==7) { $state=1 }
  print OUT "$state\n";
}
close(IN);
close(OUT);




