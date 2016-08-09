#!/usr/bin/perl
use strict;
use ProgramName;
use Fastb;

my $name=ProgramName::get();
die "$name <indir>\n" unless @ARGV==1;
my ($dir)=@ARGV;

my (%counts,%centrals);
my @files=`ls $dir`;
foreach my $file (@files) {
  chomp $file; next unless $file=~/\.fastb$/;
  my $fastb=new Fastb("$dir/$file");
  my $stateTrack=$fastb->getTrackByName("state");
  my ($begin,$end)=findPeak($stateTrack);
  my $center=int(($begin+$end)/2);
  #my $slice=$fastb;
  my $slice=$fastb->slice($begin,$end);
  my @hits;
  findHits($slice,"CTCF",\@hits);
  findHits($slice,"FOX",\@hits);
  findHits($slice,"CEBP",\@hits);
  findHits($slice,"AP1",\@hits);
  findHits($slice,"KLF",\@hits);
  findHits($slice,"GR/AR/MR",\@hits);
  foreach my $hit (@hits) {
    $counts{$hit->{factor}}++;
  }
  my $central=findCentral(\@hits,$center);
  if($central) { ++$centrals{$central->{factor} } }
}
my @factors=keys %counts;
foreach my $factor (@factors) {
  #my $count=$counts{$factor};
  my $count=$centrals{$factor};
  print "$factor\t$count\n";
}


################################################################3
sub findCentral {
  my ($hits,$center)=@_;
  my ($best,$bestDist);
  foreach my $hit (@$hits) {
    my $distance=0;
    if($hit->{begin}>$center) { $distance=$hit->{begin}-$center }
    elsif($hit->{end}<$center) { $distance=$center-$hit->{end} }
    if(!defined($bestDist) || $distance<$bestDist)
      { $best=$hit; $bestDist=$distance }
  }
  return $best;
}

sub findHits {
  my ($fastb,$factor,$hits)=@_;
  my $track=$fastb->getTrackByName($factor); die unless $track;
  my $data=$track->getData();
  my $L=@$data;
  my $begin;
  for(my $i=0 ; $i<$L ; ++$i) {
    if($i==0 && $data->[$i]>0) { $begin=$i; next }
    if($i>0 && $data->[$i]>0 && $data->[$i-1]==0) { $begin=$i }
    if($i>0 && $data->[$i]==0 && $data->[$i-1]>0) {
      push @$hits,{factor=>$factor,begin=>$begin,end=>$i};
    }
  }
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



