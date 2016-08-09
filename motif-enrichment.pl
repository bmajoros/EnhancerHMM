#!/usr/bin/perl
use strict;
use ProgramName;
use Fastb;

my $name=ProgramName::get();
die "$name <indir>\n" unless @ARGV==1;
my ($dir)=@ARGV;

my @files=`ls $dir`;
foreach my $file (@files) {
  chomp $file; next unless $file=~/\.fastb$/;
  my $fastb=new Fastb($file);
  my $stateTrack=$fastb->getTrackByName("state");
  my ($begin,$end)=findPeak($stateTrack);
  print "$begin - $end\n";
}


################################################################3
sub findPeak {
  my ($track)=@_;
  my $array=$track->getData();
  my $L=@$array;
  my ($begin,$end);
  for(my $i=0 ; $i<$L ; ++$i) { if($array->[$i]==3) {$begin=$i;last} }
  for(my $i=$begin ; $i<$L ; ++$i) { if($array->[$i]!=3) {$end=$i;last} }
  return ($begin,$end);
}



