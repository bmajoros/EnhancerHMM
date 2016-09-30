#!/usr/bin/perl
use strict;

my $BASE="/home/bmajoros/GGR/p300";
my $LS="$BASE/ls-raw.txt";


#================================================================
sub listToHash {
  my ($list)=@_;
  my $hash={};
  my $n=@$list;
  for(my $i=0 ; $i<$n ; ++$i) { $hash->{$list->[$i]}=1 }
  return $hash;
}
#================================================================
sub copyFiles {
  my ($files,$to,$notThese)=@_;
  my $n=@$files;
  my $filesCopied=0;
  for(my $i=0 ; $i<$n ; ++$i) {
    my $file=$files->[$i];
    if($notThese->{$file}) { next }
    System("cp $file $to");
    ++$filesCopied;
    last unless $filesCopied<$MAX_FILES;
  }
}
#================================================================
sub System {
  my ($cmd)=@_;
  system($cmd);
  #print "$cmd\n";
}
#================================================================

