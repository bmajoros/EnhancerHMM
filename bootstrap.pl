#!/bin/env perl
use Fastb;
$|=1;

# Globals
my $VERBOSE=0;
my $MIN_PEAK_LEN=200;
my $MAX_OFF_CENTER=500;
my $MIN_HUMP_LEN=200;
my $MIN_PEAK_DNASE=1;
my $MAX_PEAK_MARKS=0;
my $PEAK_MEAN_WIDTH=50;
my $MIN_HUMP_MARKS=0;
my $MAX_HUMP_DNASE=0;
my $MAX_FLANK=0;
my $MUMMIE=$ENV{"MUMMIE"};
my $LEFT_FLANK_STATE=1;
my $LEFT_HUMP_STATE=2;
my $PEAK_STATE=3;
my $RIGHT_HUMP_STATE=4;
my $RIGHT_FLANK_STATE=5;

# Process command line
die "bootstrap.pl <hmm> <likelihoods.txt> <outfile>\n" unless @ARGV==3;
my ($HMM,$likelihoodFile,$outfile)=@ARGV;

# Load list of files and likelihoods
my @files;
open(IN,$likelihoodFile) || die "can't open $likelihoodFile\n";
while(<IN>) {
  chomp; my @fields=split; next unless @fields>=2;
  my ($LL,$file)=@fields;
  push @files,[$LL,$file];
}
close(IN);

# Sort files by likelihood
@files=sort {$b->[0] <=> $a->[0]} @files;

# Process each file
open(OUT,">$outfile") || die "can't write to file: $outfile\n";
my $numFiles=@files;
for(my $i=0 ; $i<$numFiles ; ++$i) {
  my $file=$files[$i]->[1];
  next unless -e $file;
  next if $file=~/\.t00\./;

  # Parse and apply filters
  my $fastb=new Fastb($file);
  my $path=parse($file);
  if(!goodPath($path,$fastb)) {
    #system("/home/bmajoros/GGR/src/annotate-fastb.pl $file $HMM temp.fastb");
    #system("cat temp.fastb | $MUMMIE/fastb-to-xgraph.pl");
    next;
  }

  # Plot on screen
  #system("/home/bmajoros/GGR/src/annotate-fastb.pl $file $HMM temp.fastb");
  #system("cat temp.fastb | $MUMMIE/fastb-to-xgraph.pl");
  print OUT "$file\n";

  # Clean up
  undef $fastb;
}
close(OUT);


#=================================================================
sub parse {
  my ($file)=@_;
  my $path=[];
  open(IN,"$MUMMIE/parse $HMM $file |") || die;
  <IN>; # skip header line
  while(<IN>) {
    chomp;
    my $state=$_+0;
    push @$path,$state;
  }
  return $path;
}
#=================================================================
sub goodPath {
  my ($path,$fastb)=@_;
  my $L=@$path;

  # Parse into intervals
  my $intervals=getIntervals($path);
  my $numIntervals=@$intervals;
  die $numIntervals unless $numIntervals==5;
  my ($leftFlank,$leftHump,$peak,$rightHump,$rightFlank)=@$intervals;

  #print "$leftFlank->[0]\-$leftFlank->[1]  $leftHump->[0]\-$leftHump->[1]  $peak->[0]\-$peak->[1]  $rightHump->[0]\-$rightHump->[1]  $rightFlank->[0]\-$rightFlank->[1]\n";

  # Examine positions and widths of peak and humps
  my $peakCenter=int($peak->[0]+$peak->[1])/2;
  my $offCenter=abs($peakCenter-$L/2);
  my $peakLen=$peak->[1]-$peak->[0];
  my $leftHumpLen=$leftHump->[1]-$leftHump->[0];
  my $rightHumpLen=$rightHump->[1]-$rightHump->[0];
  if($offCenter > $MAX_OFF_CENTER ||
     $peakLen < $MIN_PEAK_LEN ||
     $leftHumpLen < $MIN_HUMP_LEN ||
     $rightHumpLen < $MIN_HUMP_LEN) {
    if($VERBOSE) { print "offCenter=$offCenter peakLen=$peakLen leftHump=$leftHumpLen rightHump=$rightHumpLen\n" }
    return 0 }

  # Compute mean values of tracks in peak
  my $peakDnase=trackIntervalMean($fastb,"DNase",$peak,$PEAK_MEAN_WIDTH);
  my $peakP300=trackIntervalMean($fastb,"EP300",$peak,$PEAK_MEAN_WIDTH);
  my $peakH3K4me1=trackIntervalMean($fastb,"H3K4me1",$peak,$PEAK_MEAN_WIDTH);
  my $peakH3K4me2=trackIntervalMean($fastb,"H3K4me2",$peak,$PEAK_MEAN_WIDTH);
  my $peakH3K27ac=trackIntervalMean($fastb,"H3K27ac",$peak,$PEAK_MEAN_WIDTH);

  # Compute mean values of tracks in left hump
  my $leftHumpDnase=trackIntervalMean($fastb,"DNase",$leftHump);
  my $leftHumpP300=trackIntervalMean($fastb,"EP300",$leftHump);
  my $leftHumpH3K4me1=trackIntervalMean($fastb,"H3K4me1",$leftHump);
  my $leftHumpH3K4me2=trackIntervalMean($fastb,"H3K4me2",$leftHump);
  my $leftHumpH3K27ac=trackIntervalMean($fastb,"H3K27ac",$leftHump);

  # Compute mean values of tracks in right hump
  my $rightHumpDnase=trackIntervalMean($fastb,"DNase",$rightHump);
  my $rightHumpP300=trackIntervalMean($fastb,"EP300",$rightHump);
  my $rightHumpH3K4me1=trackIntervalMean($fastb,"H3K4me1",$rightHump);
  my $rightHumpH3K4me2=trackIntervalMean($fastb,"H3K4me2",$rightHump);
  my $rightHumpH3K27ac=trackIntervalMean($fastb,"H3K27ac",$rightHump);

  # Compute mean values of tracks in left flank
  my $leftFlankDnase=trackIntervalMean($fastb,"DNase",$leftFlank);
  my $leftFlankP300=trackIntervalMean($fastb,"EP300",$leftFlank);
  my $leftFlankH3K4me1=trackIntervalMean($fastb,"H3K4me1",$leftFlank);
  my $leftFlankH3K4me2=trackIntervalMean($fastb,"H3K4me2",$leftFlank);
  my $leftFlankH3K27ac=trackIntervalMean($fastb,"H3K27ac",$leftFlank);

  # Compute mean values of tracks in right flank
  my $rightFlankDnase=trackIntervalMean($fastb,"DNase",$rightFlank);
  my $rightFlankP300=trackIntervalMean($fastb,"EP300",$rightFlank);
  my $rightFlankH3K4me1=trackIntervalMean($fastb,"H3K4me1",$rightFlank);
  my $rightFlankH3K4me2=trackIntervalMean($fastb,"H3K4me2",$rightFlank);
  my $rightFlankH3K27ac=trackIntervalMean($fastb,"H3K27ac",$rightFlank);

  # Apply filters
  if($peakH3K4me1>$MAX_PEAK_MARKS ||
     $peakH3K4me2>$leftHumpH3K4me2 || $peakH3K4me2>$rightHumpH3K4me2
     || $peakH3K27ac>$MAX_PEAK_MARKS || $peakDnase<$MIN_PEAK_DNASE ||
     $peakP300<$MIN_PEAK_DNASE) {
    if($VERBOSE0) { print "PEAK: $peakH3K4me1>$MAX_PEAK_MARKS || $peakH3K4me2>$MAX_PEAK_MARKS || $peakH3K27ac>$MAX_PEAK_MARKS || $peakDnase<$MIN_PEAK_DNASE || $peakP300<$MIN_PEAK_DNASE\n" }
    return 0;  }
  if($humpH3K4me1<$MIN_HUMP_MARKS || $humpH3K4me2<$MIN_HUMP_MARKS
     || $humpH3K27ac<$MIN_HUMP_MARKS || $humpDnase>$MAX_HUMP_DNASE ||
     $humpP300>$MAX_HUMP_DNASE) {
    if($VERBOSE) { print "LEFT HUMP: $humpH3K4me1<$MIN_HUMP_MARKS || $humpH3K4me2<$MIN_HUMP_MARKS || $humpH3K27ac<$MIN_HUMP_MARKS || $humpDnase>$MAX_HUMP_DNASE || $humpP300>$MAX_HUMP_DNASE\n" }
    return 0;  }
  if($humpH3K4me1<$MIN_HUMP_MARKS || $humpH3K4me2<$MIN_HUMP_MARKS
     || $humpH3K27ac<$MIN_HUMP_MARKS || $humpDnase>$MAX_HUMP_DNASE ||
     $humpP300>$MAX_HUMP_DNASE) {
    if($VERBOSE) { print "RIGHT HUMP: $humpH3K4me1<$MIN_HUMP_MARKS || $humpH3K4me2<$MIN_HUMP_MARKS || $humpH3K27ac<$MIN_HUMP_MARKS || $humpDnase>$MAX_HUMP_DNASE || $humpP300>$MAX_HUMP_DNASE\n" }
    return 0;  }
  if($leftFlankDnase>0 || $leftFlankP300>0 ||
     $leftFlankH3K4me1>$leftHumpH3K4me1 ||
     $leftFlankH3K4me2>$leftHumpH3K4me2 ||
     $leftFlankH3K27ac>$leftHumpH3K27ac) {
    if($VERBOSE) { print "LEFT FLANK: $leftFlankDnase>0 || $leftFlankP300>0 || $leftFlankH3K4me1>$leftHumpH3K4me1 || $leftFlankH3K4me2>$leftHumpH3K4me2 || $leftFlankH3K27ac>$leftHumpH3K27ac\n" }
    return 0;  }
  if($rightFlankDnase>0 || $rightFlankP300>0 ||
     $rightFlankH3K4me1>$rightHumpH3K4me1 ||
     $rightFlankH3K4me2>$rightHumpH3K4me2 ||
     $rightFlankH3K27ac>$rightHumpH3K27ac) {
    if($VERBOSE) { print "RIGHT FLANK: $rightFlankDnase>0 || $rightFlankP300>0 || $rightFlankH3K4me1>$rightHumpH3K4me1 || $rightFlankH3K4me2>$rightHumpH3K4me2 || $rightFlankH3K27ac>$rightHumpH3K27ac\n" }
    return 0;  }

  if($VERBOSE) { print "OK\n" }
  return 1;
}
#=================================================================
sub getIntervals {
  my ($path)=@_;
  my $intervals=[];
  my $L=@$path;
  my ($prev,$begin);
  for(my $i=0 ; $i<$L ; ++$i) {
    my $x=$path->[$i];
    if($x ne $prev) {
      if(defined($begin)) { push @$intervals,[$begin,$i] }
      $begin=$i;
    }
    $prev=$x;
  }
  push @$intervals,[$begin,$L];
  return $intervals;
}
#=================================================================
sub trackIntervalMean {
  my ($fastb,$trackName,$interval,$maxLen)=@_;
  my $track=$fastb->getTrackByName($trackName);
  $track=$track->slice($interval->[0],$interval->[1]);
  my $trackData=$track->getData();
  my $L=@$trackData;
  if($maxLen>0 && $L>$maxLen) {
    my $trim=$L-$maxLen;
    $track=$track->slice(int($trim/2),$L-int($trim/2));
    $trackData=$track->getData();
    $L=@$trackData;
  }
  my $sum=0;
  foreach my $x (@$trackData) { $sum+=$x }
  my $mean=$sum/$L;
  return $mean;
}
#=================================================================
#=================================================================
#=================================================================
#=================================================================
#=================================================================



