#!/bin/env perl
use strict;
use SummaryStats;
$|=1;

# Globals
my $MAX_GENES=1178;
my $BASE="/home/bmajoros/GGR/p300";
my $ACTIVATING="$BASE/increasing-new-55.txt";
my $DEACTIVATING="$BASE/decreasing-new-55.txt";
#my $ACTIVATING="$BASE/increasing-pos5-neg1.txt";
#my $DEACTIVATING="$BASE/decreasing-pos5-neg1.txt";
#my $ACTIVATING="$BASE/on-pos5-neg5.txt";
#my $DEACTIVATING="$BASE/off-pos5-neg5.txt";
my $PEAKS="$BASE/peaks.bed";
my $UP_GENES="$BASE/genes/DEGs_upreg.FDR_0.01.TSS.protein_coding.bed";
my $DOWN_GENES="$BASE/genes/DEGs_downreg.FDR_0.1.TSS.protein_coding.bed";
my $NON_DEGS="$BASE/genes/non_DEGs.FDR_0.2.TSS.protein_coding.bed";

# Load peak coordinates
my %peaks;
open(IN,$PEAKS) || die $PEAKS;
while(<IN>) {
  chomp; my @fields=split; next unless @fields>=4;
  my ($chr,$begin,$end,$id)=@fields;
  $peaks{$id} = { chr=>$chr, begin=>$begin, end=>$end } }
close(IN);
#my $numPeaks=keys %peaks;
#print "$numPeaks peaks\n";

# Load DEGs
my %DEGs; my $numGenesLoaded;
#loadGenes($UP_GENES,1,\%DEGs);
loadGenes($DOWN_GENES,-1,\%DEGs);
#loadGenes($NON_DEGS,0,\%DEGs);
#print "$numGenesLoaded genes loaded\n";

# Sort DEGs
my @chroms=keys %DEGs;
foreach my $chr (@chroms) {
  @{$DEGs{$chr}}=sort {$a->{pos} <=> $b->{pos}} @{$DEGs{$chr}}}

#my $numChr1=@{$DEGs{"chr1"}}; print "$numChr1 genes on chr1\n";

# Process activating/deactivating and compare mean distances to DEGs
my ($posDist,$posSD,$posN)=analyze($ACTIVATING,"activating-distances.txt");
my ($negDist,$negSD,$negN)=analyze($DEACTIVATING,"deactivating-distances.txt");
print "activating:   $posDist +/- $posSD N=$posN\n";
print "inactivating: $negDist +/- $negSD N=$negN\n";




#======================================================================
sub loadGenes {
  my ($GENES,$label,$hash)=@_;
  my @array;
  open(IN,$GENES) || die $GENES;
  while(<IN>) {
    chomp; my @fields=split; next unless @fields>=4;
    my ($chr,$begin,$end,$id)=@fields;
    push @array,[$chr,$begin];
  }
  close(IN);
  subsample(\@array,$MAX_GENES);
  $numGenesLoaded=@array;
  for(my $i=0 ; $i<$numGenesLoaded ; ++$i) {
    my $pair=$array[$i];
    my ($chr,$begin)=@$pair;
    push @{$hash->{$chr}},{pos=>$begin,dir=>$label};
  }
}
#======================================================================
sub analyze {
  my ($peaksFile,$outfile)=@_;
  my @distances;
  open(OUT,">$outfile") || die $outfile;
  open(IN,$peaksFile) || die $peaksFile;
  while(<IN>) {
    chomp; next unless(/\S/);
    my $record=$peaks{$_}; die $_ unless $record;
    my $chr=$record->{chr}; my $begin=$record->{begin}; my $end=$record->{end};
    next if $chr=~/_/;
    my $pos=int(($begin+$end)/2);
    next if $chr eq "chrX" || $chr eq "chrY";
    my $dist=findNearest($chr,$pos);
    if($dist>=0) {
      push @distances,$dist;
      print OUT "$dist\n";
    }
    #print "$pos\t$dist\n";
  }
  close(IN); close(OUT);
  my ($mean,$stddev,$min,$max)=SummaryStats::roundedSummaryStats(\@distances);
  my $n=@distances;
  return ($mean,$stddev,$n);
}
#======================================================================
sub shuffle {
  my ($array)=@_;
  my $n=@$array;
  for(my $i=0 ; $i<$n ; ++$i) {
    my $j=int(rand($n));
    my $tmp=$array->[$i];
    $array->[$i]=$array->[$j];
    $array->[$j]=$tmp;
  }
}
#======================================================================
sub subsample {
  my ($array,$n)=@_;
  shuffle($array);
  #my $num=@$array;  print "num=$num\n";
  splice(@$array,$n);
  #my $num=@$array;  print "num=$num\n";
}
#======================================================================
sub findNearest {
  my ($chr,$pos)=@_;
  my $A=$DEGs{$chr};
  if(!$A) { return undef } #die $chr unless $A;
  my $n=@$A;
  my $begin=0; my $end=$n-1;
  while($begin<$end) {
    my $mid=int(($begin+$end)/2);
    my $midPos=$A->[$mid]->{pos};
    if($pos<$midPos) { $end=$mid-1 }
    elsif($pos>$midPos) { $begin=$mid+1 }
    else { $begin=$end=$mid }
  }
  die unless $begin>=0 && $begin<$n;
  my $thisPos=$A->[$begin]->{pos};
  my $dist=abs($pos-$thisPos);
  if($begin>0) {
    my $prevPos=$A->[$begin-1]->{pos};
    my $prevDist=abs($pos-$prevPos);
    if($prevDist<$dist) { $dist=$prevDist }
  }
  if($begin+1<$n) {
    my $nextPos=$A->[$begin+1]->{pos};
    my $nextDist=abs($pos-$nextPos);
    if($nextDist<$dist) { $dist=$nextDist }
  }
#print "$pos $begin $end $dist\n";
  ###
  my $DEBUG=0;
  if($DEBUG) {
    my $prevPos=$A->[$begin-1]; my $thisPos=$A->[$begin];
    my $nextPos=$A->[$begin+1]; my $twoBackPos=$A->[$begin-2];
    my $twoForwardPos=$A->[$begin+2];
    my $prevDist=abs($pos-$prevPos);  my $thisDist=abs($pos-$thisPos);
    my $nextDist=abs($pos-$nextPos); my $twoBackDist=abs($pos-$twoBackPos);
    my $twoForwardDist=abs($pos-$twoForwardPos);
    if($twoBackDist<$dist || $prevDist<$dist || $thisDist<$dist ||
       $nextDist<$dist || $twoForwardDist<$dist) {
      $prevDist=log($prevDist); $twoBackDist=log($twoBackDist);
      $thisDist=log($thisDist); $twoForwardDist=log($twoForwardDist);
      $nextDist=log($nextDist); my $logDist=log($dist);
      #print "$logDist      $twoBackDist $prevDist $thisDist $nextDist $twoForwardDist\n";
      die if $prevDist<$thisDist || $nextDist<$thisDist;
    }
  }
  ###

  return $dist;
}
#======================================================================
#======================================================================
#======================================================================
#======================================================================
#======================================================================







