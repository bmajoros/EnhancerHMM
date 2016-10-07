#!/usr/bin/perl
use strict;
use ProgramName;
$|=1;

my $name=ProgramName::get();
die "$name <pos.hmm> <neg.hmm> <peaks.txt> <fastb-dir>\n" unless @ARGV==4;
my ($ACTIVE_HMM,$INACTIVE_HMM,$PEAKS,$FASTB)=@ARGV;

my $BASE="/home/bmajoros/GGR/p300";
my $MUMMIE=$ENV{"MUMMIE"};
my $PRIOR1=log(0.5);
my $PRIOR0=log(0.5);

my @times=("t00","t05","t1","t2","t3","t4","t5","t6","t7","t8","t10","t12");

open(IN,$PEAKS) || die "can't open $PEAKS\n";
while(<IN>) {
  chomp; next unless(/iter0/);
  my $peak=$_;
  my @files;
  foreach my $time (@times) {
    my $filename="$FASTB/$peak.standardized_across_all_timepoints.$time.fastb";
    next unless -e $filename;
    push @files,$filename;
  }
  next unless @files==12;
  print "$peak\t";
  foreach my $file (@files) {
    my $LL1=0+`$MUMMIE/get-likelihood $ACTIVE_HMM $file`;
    my $LL0=0+`$MUMMIE/get-likelihood $INACTIVE_HMM $file`;
    my $joint0=$LL0+$PRIOR0; my $joint1=$LL1+$PRIOR1;
    my $denom=sumLogProbs($joint0,$joint1);
    my $posterior=exp($joint1-$denom);
    $posterior=int($posterior*100+5/9)/100;
    #print "\nXXX $LL0\t$LL1\t$joint0\t$joint1\t$denom\t$posterior\n";
    print "$posterior\t"
  }
  print "\n";
}
close(IN);

sub sumLogProbs {
  my ($a,$b)=@_;
  if($a<$b) { ($a,$b)=($b,$a) }
  return $a+log(1+exp($b-$a))
}

