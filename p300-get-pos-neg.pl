#!/usr/bin/perl
use strict;

my $BASE="/home/bmajoros/GGR/p300";
my $LISTS="$BASE/active-inactive-lists";

my $activeFiles=getFiles("$LISTS/EP300_enhancers.active*txt");
my $inactiveFiles=getFiles("$LISTS/EP300_enhancers.inactive*txt");

sub getFiles {
  my ($pattern)=@_;
  my $list=[];
  my @listFiles=`ls $pattern`;
  foreach my $listFile (@listFiles) {
    chomp $listFile; next unless $listFile=~/(t\d+)/;
    my $time=$1;
    open(IN,"$LISTS/$listFile") || die "LISTS/$listFile";
    while(<IN>) {
      chomp; next unless(/iter0_peak(\d+)/);
      my $peak=$1;
      my $fastb=
	"raw/iter0_peak$peak.standardized_across_all_timepoints.$time.fastb";
      push @$list,$raw;
    }
  }
  return $list;
}


