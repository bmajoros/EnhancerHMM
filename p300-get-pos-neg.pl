#!/usr/bin/perl
use strict;

my $MAX_FILES=2000;
my $BASE="/home/bmajoros/GGR/p300";
my $LISTS="$BASE/active-inactive-lists";
my $LS="$BASE/ls-raw.txt";

my %exists;
getExists();

my $activeFiles=getFiles("$LISTS/EP300_enhancers.active*txt");
my $inactiveFiles=getFiles("$LISTS/EP300_enhancers.inactive*txt");
my $activeHash=listToHash($activeFiles);
my $inactiveHash=listToHash($inactiveFiles);
shuffle($activeFiles);
shuffle($inactiveFiles);
copyFiles($activeFiles,"pos",$inactiveHash);
copyFiles($inactiveFiles,"neg",$activeHash);

#================================================================
sub listToHash {
  my ($list)=@_;
  my $hash={};
  my $n=@$list;
  for(my $i=0 ; $i<$n ; ++$i) { $hash->{$list->[$i]}=1 }
  return $hash;
}
#================================================================
sub getExists {
  open(IN,$LS) || die $LS;
  while(<IN>) {
    chomp;
    $exists{"raw/$_"}=1;
  }
  close(IN);
}
#================================================================
sub swap {
  my ($array,$i,$j)=@_;
  my $temp=$array->[$i];
  $array->[$i]=$array->[$j];
  $array->[$j]=$temp;
}
#================================================================
sub shuffle {
  my ($list)=@_;
  my $n=@$list;
  for(my $i=0 ; $i<$n ; ++$i) { swap($list,$i,int(rand($n))) }
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
sub getFiles {
  my ($pattern)=@_;
  my $list=[];
  my @listFiles=`ls $pattern`;
  foreach my $listFile (@listFiles) {
    chomp $listFile; next unless $listFile=~/(t\d+)/;
    my $time=$1;
    open(IN,"$listFile") || die "$listFile";
    while(<IN>) {
      chomp; next unless(/iter0_peak(\d+)/);
      my $peak=$1;
      my $fastb=
	"raw/iter0_peak$peak.standardized_across_all_timepoints.$time.fastb";
      if(!$exists{$fastb}) { next }
      push @$list,$fastb;
    }
    close(IN);
  }
  return $list;
}


