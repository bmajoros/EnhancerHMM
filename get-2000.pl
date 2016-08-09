#!/usr/bin/perl
use strict;

my $MAX_COUNT=2000;
my $GGR="/home/bmajoros/GGR";

my %blacklist;
blacklist("$GGR/pos-train",\%blacklist);
blacklist("$GGR/neg-train",\%blacklist);

process("$GGR/nonredundant/high","$GGR/2000/high");
process("$GGR/nonredundant/bg","$GGR/2000/bg");


####################################################################

sub System {
  my ($cmd)=@_;
  #print "$cmd\n";
  system($cmd);
}

sub process {
  my ($fromDir,$toDir)=@_;
  my @files=`ls $fromDir`;
  my $numFiles=@files;
  my $count=0;
  for(my $i=0 ; $i<$numFiles ; ++$i) {
    my $file=$files[$i]; chomp $file;
    next unless $file=~/\.fastb$/;
    $file=~/peak(\d+)\.t(\d+)/ || die $file;
    my ($peak,$time)=($1,$2);
    next unless $time>6;
    my $key="$peak $time";
    next if $blacklist{$key};
    System("cp $fromDir/$file $toDir");
    ++$count;
    if($count>=$MAX_COUNT) { return }
  }
}

sub blacklist {
  my ($dir,$hash)=@_;
  my @files=`ls $dir`;
  my $numFiles=@files;
  for(my $i=0 ; $i<$numFiles ; ++$i) {
    my $file=$files[$i]; chomp $file;
    next unless $file=~/\.fastb$/;
    $file=~/peak(\d+)\.t(\d+)/ || die $file;
    my ($peak,$time)=($1,$2);
    my $key="$peak $time";
    $hash->{$key}=1;
  }
}





