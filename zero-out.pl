#!/usr/bin/perl
use strict;

my $KILL="H3K4me2";

system("m/fastb-zero-out-tracks.pl pos-test $KILL pos-temp");
system("m/fastb-zero-out-tracks.pl neg-test $KILL neg-temp");
system("src/classify.pl pos-temp neg-temp model1/trained-pos5.hmm model1/trained-neg1.hmm");
system("roc.pl tmp.roc > roc.txt");
system("area-under-ROC.pl roc.txt");




