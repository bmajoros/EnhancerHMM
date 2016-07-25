#!/usr/bin/perl
use strict;

# Global settings
my $BASE="/home/bmajoros/GGR";
my $MODEL_DIR="$BASE/model1";
my $SCHEMA="$MODEL_DIR/model1.schema";
my $MUMMIE=$ENV{"MUMMIE"};
my $OUTFILE="$MODEL_DIR/pos7.hmm";

system("cd $MODEL_DIR ; ../m/hmm-extract-state trained-pos5.hmm all state");

# Make a metamodel that links the individual states into a chain
my $metamodel="$MODEL_DIR/pos7.metamodel";
open(META,">$metamodel") || die;
print META "0 -> 1 : 1\n";
for(my $i=1 ; $i<7 ; ++$i) {
  my $next=$i+1;
  print META "$i ->$i : 0.997\n";
  print META "$i -> $next : 0.003\n" }
print META "7 -> 7 : 0.997\n";
print META "7 -> 0 : 0.003\n";
close(META);

# Make submodel file that lists the state filenames
my $submodels="$MODEL_DIR/pos7.submodels";
open(SUB,">$submodels") || die;
print SUB "1 = $MODEL_DIR/state1.hmm\n";
print SUB "2 = $MODEL_DIR/state2.hmm\n";
print SUB "3 = $MODEL_DIR/state3.hmm\n";
print SUB "4 = $MODEL_DIR/state4.hmm\n";
print SUB "5 = $MODEL_DIR/state3.hmm\n";
print SUB "6 = $MODEL_DIR/state4.hmm\n";
print SUB "7 = $MODEL_DIR/state5.hmm\n";
close(SUB);

# Combine states into one model
system("$MUMMIE/model-combiner $metamodel $submodels $OUTFILE");

# Clean up
for(my $i=1 ; $i<=5 ; ++$i) { system("rm $MODEL_DIR/state$i.hmm") }




