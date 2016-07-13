#!/usr/bin/perl
use strict;
use ProgramName;

my $name=ProgramName::get();
die "$name <indir> <outdir-train> <outdir-test>\n" unless @ARGV==3;
my ($indir,$outTrain,$outTest)=@ARGV;





