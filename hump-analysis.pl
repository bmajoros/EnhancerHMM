#!/bin/env perl
use strict;
use ProgramName;

my $TOLERANCE=20;
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
  my @signs;
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
    if($diff>$TOLERANCE) { $sign='>'}
    elsif($diff<-$TOLERANCE) { $sign='<' }
    push @signs,$sign;
  }
  my ($flips,$randomFlips)=getFlips(\@signs);
  print "$flips\t$randomFlips\n";
}
close(TRAJ);



#=================================================================
sub randomize {
  my ($array)=@_;
  my $random=[];
  my $n=@$array;
  for(my $i=0 ; $i<$n ; ++$i) { $random->[$i]=rand(1)<0.5 ? '<' : '>'}
  return $random;
}
#=================================================================
sub shuffle {
  my ($array)=@_;
  my $random=[];
  @$random=@$array;
  my $n=@$random;
  for(my $i=0 ; $i<$n ; ++$i) {
    my $j=int(rand($n));
    my $tmp=$array->[$i];
    $array->[$i]=$array->[$j];
    $array->[$j]=$tmp;
  }
  return $random;
}
#=================================================================
sub getFlips {
  my ($signs)=@_;
  #my $random=shuffle($signs);
  my $random=randomize($signs);
  my $flips=countFlips($signs);
  my $randomFlips=countFlips($random);
  return ($flips,$randomFlips);
}
#=================================================================
sub countFlips {
  my ($signs)=@_;
  my $n=@$signs;
  my $prev; my $flips=0;
  for(my $i=0 ; $i<$n ; ++$i) {
    my $sign=$signs->[$i];
    if($sign ne "=") {
      if(defined($prev)) { if($sign ne $prev) {++$flips} }
      $prev=$sign;
    }
  }
  return $flips;
}
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
