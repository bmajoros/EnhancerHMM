#!/usr/bin/perl
use strict;

# Global settings
my $NUM_COMPONENTS=3;
my $BASE="/home/bmajoros/GGR";
my $MODEL_DIR="$BASE/model1";
my $SCHEMA="$MODEL_DIR/model1.schema";
my $MUMMIE=$ENV{"MUMMIE"};

# Make the TGF file
system("$MUMMIE/make-tgf.pl $SCHEMA full $MODEL_DIR/tgf.tgf");

# Make the linear HMM
make5stateHMM("$MODEL_DIR/pos5.hmm");

# Make a single-state HMM as baseline


##################################################################
sub make5stateHMM
{
  my ($outfile)=@_;

  # Make individual states
  my $state1=makeLoopingState("$MODEL_DIR/fg1.hmm");
  my $state2=makeLoopingState("$MODEL_DIR/fg2.hmm");
  my $state3=makeLoopingState("$MODEL_DIR/fg3.hmm");
  my $state4=makeLoopingState("$MODEL_DIR/fg4.hmm");
  my $state5=makeLoopingState("$MODEL_DIR/fg5.hmm");

  # Make a metamodel that links the individual states into a chain
  my $metamodel="$MODEL_DIR/pos5.metamodel";
  open(META,">$metamodel") || die;
  print META "0 -> 1 : 1\n";
  for(my $i=1 ; $i<5 ; ++$i) {
    my $next=$i+1;
    print META "$i -> $next : 1\n";
  }
  print META "5 -> 0 : 1\n";
  close(META);

  # Make submodel file that lists the state filenames
  my $submodels="$MODEL_DIR/pos5.submodels";
  open(SUB,">$submodels") || die;
  for(my $i=1 ; $i<=5 ; ++$i) { print SUB "$i = $MODEL_DIR/fg$i.hmm\n" }
  close(SUB);

  # Combine states into one model
  system("$MUMMIE/model-combiner $metamodel $submodels $outfile");

  # Set emissions to a good initial value to guide EM in the right direction
  system("$MUMMIE/hmm-edit $outfile MIX 1 0 1 MIX 2 1 1 MIX 3 2 1 MIX 4 1 1 MIX 5 0 1");
  system("");
}



sub makeLoopingState
{
  my ($filename)=@_;
  my $topology="0";
  my $order=0;
  my $connectedness=1;
  system("$MUMMIE/random-HMM -c $topology 2 $connectedness $NUM_COMPONENTS $SCHEMA $order $filename");
}


