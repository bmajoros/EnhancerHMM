#!/bin/env perl
use strict;
use ProgramName;

my $BASE="/home/bmajoros/GGR/p300";
my $FASTB="$BASE/no-dna";
my $MUMMIE=$ENV{"MUMMIE"};
my @times=("t00","t05","t1","t2","t3","t4","t5","t6","t7","t8","t10","t12");

my $name=ProgramName::get();
die "$name <trajectories.txt> <pos.hmm>\n" unless @ARGV==2;
my ($trajectoriesFile,$hmmFile)=@ARGV;

open(TRAJ,$trajectoriesFile) || die "can't open $trajectoriesFile\n";
while(<TRAJ>) {
  chomp; my @fields=split; next unless @fields>=13;
  my $peak=shift @fields;
  my @ones;
  my $numFields=@fields; die unless $numFields==12;
  for(my $i=0 ; $i<$numFields ; ++$i) {
    my $posterior=$fields[$i];
    if($posterior==1) { push @ones,$i }
  }
  next unless @ones>=12;
  print "$peak\t";
  foreach my $i (@ones) {
    my $time=$times[$i];
    my $fastb="$FASTB/$peak.standardized_across_all_timepoints.$time.fastb";
    my $path=parse($fastb,$hmmFile);
    my $intervals=getIntervals($path);
    die unless @$intervals==5;
    my ($leftFlank,$leftHump,$peak,$rightHump,$rightFlank)=@$intervals;
    my $leftHumpLen=$leftHump->[1]-$leftHump->[0];
    my $rightHumpLen=$rightHump->[1]-$rightHump->[0];
    my $diff=$rightHumpLen-$leftHumpLen;
    my $sign='=';
    if($diff>20) { $sign='>'}
    elsif($diff<-20) { $sign='<' }
    print "$sign  ";
  }
  print "\n";
}
close(TRAJ);



#=================================================================
sub parse {
  my ($file,$HMM)=@_;
  my $path=[];
  open(PARSE,"$MUMMIE/parse $HMM $file |") || die;
  <PARSE>; # skip header line
  while(<PARSE>) {
    chomp;
    my $state=$_+0;
    push @$path,$state;
  }
  return $path;
}
#=================================================================
sub getIntervals {
  my ($path)=@_;
  my $intervals=[];
  my $L=@$path;
  my ($prev,$begin);
  for(my $i=0 ; $i<$L ; ++$i) {
    my $x=$path->[$i];
    if($x ne $prev) {
      if(defined($begin)) { push @$intervals,[$begin,$i] }
      $begin=$i;
    }
    $prev=$x;
  }
  push @$intervals,[$begin,$L];
  return $intervals;
}
