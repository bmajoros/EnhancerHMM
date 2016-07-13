#!/usr/bin/perl
use strict;

my $NUM_LINEAR_STATES=5;
my $NUM_COMPONENTS=3;
my $BASE="/home/GGR";
my $MODEL_DIR="$BASE/model1";
my $SCHEMA="$MODEL_DIR/model1.schema";

sub makeLinearHMM
{
  my $topology="0";
  my $order=0;
  my $numStates=2;
  my $connectedness=1;
  system("$MUMMIE/random-HMM -c $topology $numStates $connectedness $NUM_COMPONENTS $SCHEMA $order model1/state.hmm");


}


