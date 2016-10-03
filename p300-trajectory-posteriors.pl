#!/usr/bin/perl
use strict;

my $BASE="/home/bmajoros/GGR/p300";
my $MODEL="$BASE/model";
my $ACTIVE_HMM="$MODEL/trained-pos5.hmm"
my $INACTIVE_HMM="$MODEL/trained-neg1.hmm";
my $PEAKS="$BASE/peaks.txt";
my $FASTB="$BASE/no-dna";
my $PRIOR1=log(0.5);
my $PRIOR0=log(0.5);

my @times=("t00","t05","t1","t2","t3","t4","t5","t6","t7","t8","t9","t10",
	   "t11","t12");

open(IN,$PEAKS) || die "can't open $PEAKS\n";
while(<IN>) {
  chomp; next unless(/iter0/);
  my $peak=$_;
  my @files;
  foreach my $time (@times) {
    my $filename="$FASTB/$peak.standardized_across_all_timepoints.$time.fastb";
    next unless -e $filename;
    push @files,$file;
  }
  next unless @files==13;
  print "$peak\t";
  foreach my $file (@files) {
    my $LL1=0+`$MUMMIE/get-likelihood $ACTIVE_HMM $file`;
    my $LL0=0+`$MUMMIE/get-likelihood $INACTIVE_HMM $file`;
    my $denom=sumLogProbs($LL0+$PRIOR0,$LL1+$PRIOR1);
    my $posterior=$LL1+$PRIOR1/$denom;
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

