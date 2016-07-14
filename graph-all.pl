#!/usr/bin/perl
use strict;
use ProgramName;
$|=1;

my $MUMMIE=$ENV{"MUMMIE"};

my $name=ProgramName::get();
die "$name <dir>\n" unless @ARGV==1;
my ($dir)=@ARGV;

my @files=`ls $dir`;
my $n=@files;
for(my $i=0 ; $i<$n ; ++$i) {
  my $file=$files[$i]; chomp $file;
  print "$dir/$file  ---  press ENTER...\n";
  system("cat $dir/$file | $MUMMIE/fastb-to-xgraph.pl");
  <STDIN>;

}
