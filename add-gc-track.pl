#!/usr/bin/perl
use strict;
use ProgramName;
use Fastb;

my $name=ProgramName::get();
die "$name <in.fastb> <window-size> <out.fastb>\n" unless @ARGV==3;
my ($infile,$windowSize,$outfile)=@ARGV;

my $fastb=new Fastb($infile);
my $L=$fastb->getLength();
my $dnaTrack=$fastb->getTrackByName("DNA");
die "no DNA track found" unless $dnaTrack;


my $GC=[];
my $lastWindow=$L-$windowSize;
for(my $i=0 ; $i<=$lastWindow ; ++$i) {
  
}
my $gcTrack=new FastbTrack("continuous","GC",$GC);



