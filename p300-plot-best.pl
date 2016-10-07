#!/bin/env perl
$|=1;

my $BASE="/home/bmajoros/GGR/p300";
#my $LIKELIHOODS="$BASE/likelihoods.txt";
my $LIKELIHOODS="$BASE/bootstrap-likelihoods.txt";
my $RETRAIN="$BASE/retrain.txt";
my $DIR="no-dna";
my $MUMMIE=$ENV{"MUMMIE"};
my $HMM="$BASE/model/bootstrap-pos5.hmm";

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
  #next unless $retrain{$file};
  push @files,[$LL,$file];
}
close(IN);

# Sort files by likelihood
#@files=sort {$b->[0] <=> $a->[0]} @files;
@files=sort {$a->[0] <=> $b->[0]} @files;

my $numFiles=@files;
for(my $i=0 ; $i<$numFiles ; ++$i) {
  my $pair=$files[$i];
  my ($LL,$file)=@$pair;
  print "$LL\t$file\n";
  system("src/annotate-fastb.pl $DIR/$file $HMM temp.fastb");
  system("cat temp.fastb | $MUMMIE/fastb-to-xgraph.pl");
  #system("cat $DIR/$file | $MUMMIE/fastb-to-xgraph.pl");
}


