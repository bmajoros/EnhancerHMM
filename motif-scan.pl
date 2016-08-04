#!/usr/bin/perl
use strict;
use ProgramName;
use TempFilename;
use Fastb;

my $THRESHOLD=0.0001;
my $MARKOV="/data/reddylab/software/meme_4.10.0/bin/fasta-get-markov";
my $MAST="/data/reddylab/software/meme_4.10.0/bin/mast";
my $EXTRACT="/home/bmajoros/MUMMIE/fastb-to-fasta.pl";
my $TEMP_FASTA=TempFilename::generate();
my $TEMP_BACKGROUND=TempFilename::generate();
my $TEMP_OUT=TempFilename::generate();

my $name=ProgramName::get();
die "$name <in-fastb-dir> <out-fastb-dir>" unless @ARGV==2;
my ($inDir,$outDir)=@ARGV;

my @files=`ls $inDir`;
my $n=@files;
for(my $i=0 ; $i<$n ; ++$i) {
  my $file=$files[$i]; chomp $file;
  next unless $file=~/\.fastb$/;
  system("$EXTRACT $inDir/$file $TEMP_FASTA");
  system("$MARKOV -m 1 $TEMP_FASTA > $TEMP_BACKGROUND");
  my $fastb=new Fastb("$indir/$file");
  my $L=$fastb->getLength();
  my $motifCode="0"x$L;

  ###
  my $MOTIF="main-motifs/GR.meme";
  ###

  system("$MAST $MOTIF $TEMP_FASTA -hit_list -mt $THRESHOLD -bfile $TEMP_BACKGROUND > $TEMP_OUT");
  open(IN,$TEMP_OUT) || die $TEMP_OUT;
  while(<IN>) {
    next if(/^#/);
    chomp; my @fields=split; next unless @fields>=6;
    my ($name,$strand,$begin,$end,$score,$P)=@fields;
    for(my $i=$begin ; $i<$end ; ++$i) { substr($motifCode,$i,1,"1") }
  }
  close(IN);

# # All non-overlapping hits in all sequences from "tmp.fasta".
# # sequence_name motif hit_start hit_end score hit_p-value
# iter0_peak130864.t05.standardized -1 9 25 -1904.73 1.14e-02
# iter0_peak130864.t05.standardized +1 33 49 -1094.17 2.93e-03
# iter0_peak130864.t05.standardized +1 55 71 -2020.52 1.36e-02
# iter0_peak130864.t05.standardized -1 75 91 -2252.10 1.90e-02
# iter0_peak130864.t05.standardized +1 148 164 -1765.77 9.23e-03
# # mast main-motifs/GR.meme tmp.fasta -hit_list -mt 0.05 -bfile tmp.bg


  exit;
}

unlink($TEMP_FASTA); unlink($TEMP_BACKGROUND); unlink($TEMP_OUT);


