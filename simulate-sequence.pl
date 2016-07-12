#!/usr/bin/perl
use strict;
use ProgramName;
use Fastb;

my @ALPHABET=("A","C","G","T");

my $name=ProgramName::get();
die "$name <in.fastb> <out.fastb>\n" unless @ARGV==2;
my ($infile,$outfile)=@ARGV;

my $fastb=new Fastb($infile);
my $L=$fastb->getLength();
my $seq;
for(my $i=0 ; $i<$L ; ++$i) { $seq.=$ALPHABET[int(rand(4))] }
my $track=new FastbTrack("discrete","DNA",$seq,"/source=simulated");
$fastb->addTrack($track);
$fastb->save($outfile);




