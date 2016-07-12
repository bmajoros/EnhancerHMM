#!/usr/bin/perl
use strict;
use ProgramName;
use Fastb;

my $name=ProgramName::get();
die "$name <in.fastb> <window-size> <out.fastb>\n" unless @ARGV==3;
my ($infile,$windowSize,$outfile)=@ARGV;
my $halfWindow=int($windowSize/2);

my $fastb=new Fastb($infile);
my $L=$fastb->getLength();
my $dnaTrack=$fastb->getTrackByName("DNA");
my $seq=$dnaTrack->getData();
$seq="\U$seq";
die "no DNA track found" unless $dnaTrack;

my $GC=[];
my $lastWindow=$L-$windowSize;
for(my $pos=0 ; $pos<=$lastWindow ; ++$pos) {
  my $gc=getGC($pos);
  my $center=$pos+$halfWindow;
  $GC->[$center]=$gc;
  if($pos==0) { install(0,$center,$gc,$GC) }
  elsif($pos==$lastWindow) { install($center+1,$L,$gc,$GC) }
}
my $gcTrack=new FastbTrack("continuous","GC",$GC);
$fastb->addTrack($gcTrack);
$fastb->dropTrack("DNA");
$fastb->save($outfile);


sub install
{
  my ($from,$to,$value,$array)=@_;
  for(my $i=$from ; $i<$to ; ++$i) { $array->[$i]=$value }
}



sub getGC
{
  my ($pos)=@_;
  my $substr=substr($seq,$pos,$windowSize);
  my $gc=$substr=~s/([GC])/$1/g;
  my $acgt=$substr=~s/([ACGT])/$1/g;
  return $gc/$acgt;
}


