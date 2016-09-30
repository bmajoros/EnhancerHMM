#!/usr/bin/perl
use strict;

# Global settings
my $NUM_COMPONENTS=3;
my $BASE="/home/bmajoros/GGR/p300";
my $MODEL_DIR="$BASE/model";
my $SCHEMA="$MODEL_DIR/p300.schema";
my $MUMMIE=$ENV{"MUMMIE"};
my $DNase=0;
my $P300=1
my $H3K27ac=2;
my $H3K4me1=3;
my $H3K4me2=4;

# Make the TGF file
system("$MUMMIE/make-tgf.pl $SCHEMA full $MODEL_DIR/tgf.tgf");

# Make the linear HMM
make5stateHMM("$MODEL_DIR/pos5.hmm");

# Make a single-state HMM as baseline
#make1stateHMM("$MODEL_DIR/pos1.hmm");  # uses mixture model
simple1StateHMM("$MODEL_DIR/pos1.hmm"); # no mixture model

# Make a background 5-state HMM
#make5stateBackground("$MODEL_DIR/neg5.hmm");

# Make a background 1-state HMM
#make1stateHMM("$MODEL_DIR/neg1.hmm"); # uses mixture model
simple1StateHMM("$MODEL_DIR/neg1.hmm"); # no mixture model


##################################################################
##################################################################
##################################################################


sub simple1StateHMM # just 1 mixture component
{
  my ($filename)=@_;
  my $topology="0";
  my $order=0;
  my $connectedness=1;
  my $mixtureComponents=1;
  system("$MUMMIE/random-HMM -c $topology 2 $connectedness $mixtureComponents $SCHEMA $order $filename");
}



sub make1stateHMM # uses a 3-part mixture
{
  my ($outfile)=@_;
  my $state=makeLoopingState($outfile);
}



sub make5stateHMM
{
  my ($outfile)=@_;

  # Make individual states
  makeLoopingState("$MODEL_DIR/fg1.hmm");
  makeLoopingState("$MODEL_DIR/fg2.hmm");
  makeLoopingState("$MODEL_DIR/fg3.hmm");
  makeLoopingState("$MODEL_DIR/fg4.hmm");
  makeLoopingState("$MODEL_DIR/fg5.hmm");

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
  system("$MUMMIE/hmm-edit $outfile MIX 1 0 1 MIX 1 1 0 MIX 1 2 0");
  system("$MUMMIE/hmm-edit $outfile MIX 2 0 0 MIX 2 1 1 MIX 2 2 0");
  system("$MUMMIE/hmm-edit $outfile MIX 3 0 0 MIX 3 1 0 MIX 3 2 1");
  system("$MUMMIE/hmm-edit $outfile MIX 4 0 0 MIX 4 1 1 MIX 4 2 0");
  system("$MUMMIE/hmm-edit $outfile MIX 5 0 1 MIX 5 1 0 MIX 5 2 0");
  for(my $track=0 ; $track<4 ; ++$track)
    { system("$MUMMIE/hmm-edit $outfile MEAN 0 $track 0") }
  system("$MUMMIE/hmm-edit $outfile MEAN 1 $DNase 0");
  system("$MUMMIE/hmm-edit $outfile MEAN 1 $H3K27ac 1");
  system("$MUMMIE/hmm-edit $outfile MEAN 1 $H3K4me1 1");
  system("$MUMMIE/hmm-edit $outfile MEAN 1 $H3K4me2 1");
  system("$MUMMIE/hmm-edit $outfile MEAN 2 $DNase 2");
  system("$MUMMIE/hmm-edit -- $outfile MEAN 2 $H3K27ac -2");
  system("$MUMMIE/hmm-edit -- $outfile MEAN 2 $H3K4me1 -2");
  system("$MUMMIE/hmm-edit -- $outfile MEAN 2 $H3K4me2 -2");

  # Clean up
  for(my $i=1 ; $i<=5 ; ++$i) { system("rm $MODEL_DIR/fg$i.hmm") }
}



sub make5stateBackground
{
  my ($outfile)=@_;

  # Make individual states
  makeLoopingState("$MODEL_DIR/bg1.hmm");
  makeLoopingState("$MODEL_DIR/bg2.hmm");
  makeLoopingState("$MODEL_DIR/bg3.hmm");
  makeLoopingState("$MODEL_DIR/bg4.hmm");
  makeLoopingState("$MODEL_DIR/bg5.hmm");

  # Make a metamodel that links the individual states into a chain
  my $metamodel="$MODEL_DIR/neg5.metamodel";
  open(META,">$metamodel") || die;
  print META "0 -> 1 : 1\n";
  for(my $i=1 ; $i<5 ; ++$i) {
    my $next=$i+1;
    print META "$i -> $next : 1\n";
  }
  print META "5 -> 0 : 1\n";
  close(META);

  # Make submodel file that lists the state filenames
  my $submodels="$MODEL_DIR/neg5.submodels";
  open(SUB,">$submodels") || die;
  for(my $i=1 ; $i<=5 ; ++$i) { print SUB "$i = $MODEL_DIR/bg$i.hmm\n" }
  close(SUB);

  # Combine states into one model
  system("$MUMMIE/model-combiner $metamodel $submodels $outfile");

  # Clean up
  for(my $i=1 ; $i<=5 ; ++$i) { system("rm $MODEL_DIR/bg$i.hmm") }
}



sub makeLoopingState
{
  my ($filename)=@_;
  my $topology="0";
  my $order=0;
  my $connectedness=1;
  system("$MUMMIE/random-HMM -c $topology 2 $connectedness $NUM_COMPONENTS $SCHEMA $order $filename");
}


