#!/bin/env perl
use Fastb;
$|=1;

# Globals
my $MIN_PEAK_LEN=200;
my $MAX_OFF_CENTER=500;
my $MIN_HUMP_LEN=200;
my $MIN_PEAK_DNASE=1;
my $MAX_PEAK_MARKS=0;
my $MIN_HUMP_MARKS=0.5;
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
  next unless goodPath($path,$fastb);

  # Plot on screen
  system("/home/bmajoros/GGR/src/annotate-fastb.pl $file $HMM temp.fastb");
  system("cat temp.fastb | $MUMMIE/fastb-to-xgraph.pl");
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

  print "$leftFlank->[0]\-$leftFlank->[1]  $leftHump->[0]\-$leftHump->[1]  $peak->[0]\-$peak->[1]  $rightHump->[0]\-$rightHump->[1]  $rightFlank->[0]\-$rightFlank->[1]\n";

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
    print "offCenter=$offCenter peakLen=$peakLen leftHump=$leftHumpLen rightHump=$rightHumpLen\n";
    return 0 }

  # Apply filters to values of tracks in the peak
  my $peakDnase=trackIntervalMean($fastb,"DNase",$peak);
  my $peakP300=trackIntervalMean($fastb,"EP300",$peak);
  my $peakH3K4me1=trackIntervalMean($fastb,"H3K4me1",$peak);
  my $peakH3K4me2=trackIntervalMean($fastb,"H3K4me2",$peak);
  my $peakH3K27ac=trackIntervalMean($fastb,"H3K27ac",$peak);
  if($peakH3K4me1>$MAX_PEAK_MARKS || $peakH3K4me2>$MAX_PEAK_MARKS
     || $peakH3K27ac>$MAX_PEAK_MARKS || $peakDnase<$MIN_PEAK_DNASE ||
     $peakP300<$MIN_PEAK_DNASE) {
    print "PEAK: $peakH3K4me1>$MAX_PEAK_MARKS || $peakH3K4me2>$MAX_PEAK_MARKS || $peakH3K27ac>$MAX_PEAK_MARKS || $peakDnase<$MIN_PEAK_DNASE || $peakP300<$MIN_PEAK_DNASE\n";
    return 0;  }

  # Apply filters to values of tracks in the left hump
  my $humpDnase=trackIntervalMean($fastb,"DNase",$leftHump);
  my $humpP300=trackIntervalMean($fastb,"EP300",$leftHump);
  my $humpH3K4me1=trackIntervalMean($fastb,"H3K4me1",$leftHump);
  my $humpH3K4me2=trackIntervalMean($fastb,"H3K4me2",$leftHump);
  my $humpH3K27ac=trackIntervalMean($fastb,"H3K27ac",$leftHump);
  if($humpH3K4me1<$MIN_HUMP_MARKS || $humpH3K4me2<$MIN_HUMP_MARKS
     || $humpH3K27ac<$MIN_HUMP_MARKS || $humpDnase>$MAX_HUMP_DNASE ||
     $humpP300>$MAX_HUMP_DNASE) {
    print "LEFT HUMP: $humpH3K4me1<$MIN_HUMP_MARKS || $humpH3K4me2<$MIN_HUMP_MARKS || $humpH3K27ac<$MIN_HUMP_MARKS || $humpDnase>$MAX_HUMP_DNASE || $humpP300>$MAX_HUMP_DNASE\n";
    return 0;  }

  # Apply filters to values of tracks in the right hump
  my $humpDnase=trackIntervalMean($fastb,"DNase",$rightHump);
  my $humpP300=trackIntervalMean($fastb,"EP300",$rightHump);
  my $humpH3K4me1=trackIntervalMean($fastb,"H3K4me1",$rightHump);
  my $humpH3K4me2=trackIntervalMean($fastb,"H3K4me2",$rightHump);
  my $humpH3K27ac=trackIntervalMean($fastb,"H3K27ac",$rightHump);
  if($humpH3K4me1<$MIN_HUMP_MARKS || $humpH3K4me2<$MIN_HUMP_MARKS
     || $humpH3K27ac<$MIN_HUMP_MARKS || $humpDnase>$MAX_HUMP_DNASE ||
     $humpP300>$MAX_HUMP_DNASE) {
    print "RIGHT HUMP: $humpH3K4me1<$MIN_HUMP_MARKS || $humpH3K4me2<$MIN_HUMP_MARKS || $humpH3K27ac<$MIN_HUMP_MARKS || $humpDnase>$MAX_HUMP_DNASE || $humpP300>$MAX_HUMP_DNASE\n";
    return 0;  }

  print "OK\n";
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
  my ($fastb,$trackName,$interval)=@_;
  my $track=$fastb->getTrackByName($trackName);
  $track->slice($interval->[0],$interval->[1]);
  my $trackData=$track->getData();
  my $L=@$trackData;
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



