#!/usr/bin/perl
use strict;

my $GGR="/home/bmajoros/GGR";

my %blacklist;
blacklist("$GGR/pos-train",\%blacklist);
blacklist("$GGR/neg-train",\%blacklist);

process("$GGR/nonredundant/high","$GGR/2000/high");
process("$GGR/nonredundant/bg","$GGR/2000/bg");


####################################################################

sub System {
  my ($cmd)=@_;
  print "$cmd\n";
  #system($cmd);
}

sub process {
  my ($fromDir,$toDir)=@_;
  my @files=`ls $fromDir`;
  my $numFiles=@files;
  for(my $i=0 ; $i<$numFiles ; ++$i) {
    my $file=$files[$i]; chomp $file;
    next unless $file=~/\.fastb$/;
    $file=~/peak(\d+)\.t(\d+)/ || die $fastb;
    my ($peak,$time)=($1,$2);
    my $key="$peak $time";
    next if $hash->{$key};
    System("cp $fromDir/$file $toDir");
  }
}

sub blacklist {
  my ($dir,$hash)=@_;
  my @files=`ls $dir`;
  my $numFiles=@files;
  for(my $i=0 ; $i<$numFiles ; ++$i) {
    my $file=$files[$i]; chomp $file;
    next unless $file=~/\.fastb$/;
    $file=~/peak(\d+)\.t(\d+)/ || die $fastb;
    my ($peak,$time)=($1,$2);
    my $key="$peak $time";
    $hash->{$key}=1;
  }
}





