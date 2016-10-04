#!/bin/env perl
$|=1;

my $MAX_FILES=2000;
my $BASE="/home/bmajoros/GGR/p300";
my $OUTDIR="$BASE/best";
my $LIKELIHOODS="$BASE/likelihoods.txt";
my $RETRAIN="$BASE/retrain.txt";
my $MUMMIE=$ENV{"MUMMIE"};

# Load list of files for retraining
my %retrain;
open(IN,$RETRAIN) || die "can't open $RETRAIN\n";
while(<IN>) {
  chomp;
  $retrain{$_}=1;
}
close(IN);

# Load list of files and likelihoods
my @files;
open(IN,$LIKELIHOODS) || die "can't open $LIKELIHOODS\n";
while(<IN>) {
  chomp; my @fields=split; next unless @fields>=2;
  my ($LL,$file)=@fields;
  next unless $retrain{$file};
  push @files,[$LL,$file];
}
close(IN);

# Sort files by likelihood
@files=sort {$b->[0] <=> $a->[0]} @files;

my $numFiles=@files;
if($numFiles>$MAX_FILES) { $numFiles=$MAX_FILES }
for(my $i=0 ; $i<$numFiles ; ++$i) {
  my $pair=$files[$i];
  my ($LL,$file)=@$pair;
  system("cp $file $OUTDIR");
}


