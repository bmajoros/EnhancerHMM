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
  my $fastb=new Fastb("$dir/$file");
  my $stateTrack=$fastb->getTrackByName("state");
  my ($begin,$end)=findPeak($stateTrack);
  my $center=int(($begin+$end)/2);
  my $slice=$fastb->slice($begin,$end);
  my @hits;
  findHits($slice,"CTCF",\@hits);
  findHits($slice,"FOX",\@hits);
  findHits($slice,"CEBP",\@hits);
  findHits($slice,"AP1",\@hits);
  findHits($slice,"KLF",\@hits);
  findHits($slice,"GR/AR/MR",\@hits);
}


################################################################3
sub findHits {
  my ($fastb,$factor,$hits)=@_;
  my $track=$fastb->getTrackByName($factor); die unless $track;
  
}

sub findPeak {
  my ($track)=@_;
  my $array=$track->getData();
  my $L=@$array;
  my ($begin,$end);
  for(my $i=0 ; $i<$L ; ++$i) { if($array->[$i]==3) {$begin=$i;last} }
  for(my $i=$begin ; $i<$L ; ++$i) { if($array->[$i]!=3) {$end=$i;last} }
  return ($begin,$end);
}



