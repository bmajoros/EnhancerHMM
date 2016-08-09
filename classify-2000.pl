#!/usr/bin/perl
use strict;
use ProgramName;
use Fastb;

my $THRESHOLD=1000;
my $MUMMIE=$ENV{"MUMMIE"};

my $name=ProgramName::get();
die "$name <from-dir> <to-dir> <pos.hmm> <neg.hmm>\n" unless @ARGV==4;
my ($fromDir,$toDir,$posHMM,$negHMM)=@ARGV;

process($fromDir,$toDir,$posHMM,$negHMM,1);

#########################################################################
sub process
{
  my ($fromDir,$toDir,$posHMM,$negHMM,$label)=@_;
  my @files=`ls $fromDir`;
  my $n=@files;
  for(my $i=0 ; $i<$n ; ++$i) {
    my $file=$files[$i]; chomp $file;
    next unless $file=~/\.fastb$/;
    my $numer=`$MUMMIE/get-likelihood $posHMM $fromDir-noDNA/$file`;
    my $denom=`$MUMMIE/get-likelihood $negHMM $fromDir-noDNA/$file`;
    my $ratio=$numer-$denom;
    #print "$ratio\n";
    if($ratio>=$THRESHOLD) {
      System("/home/bmajoros/GGR/src/annotate-fastb.pl $fromDir-noDNA/$file $posHMM $toDir/$file");
      my $fastb=new Fastb("$fromDir/$file");
      my $track=$fastb->getTrackByName("dna");
      my $string=$track->getData();
      open(OUT,">>$toDir/$file");
      print OUT ">dna\n$string\n";
      close(OUT);
    }
  }
}

sub System {
  my ($cmd)=@_;
  #print "$cmd\n";
  system($cmd);
}








