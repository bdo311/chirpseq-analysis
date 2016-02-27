#!/usr/bin/perl

# this script is used to change .wig file to a fake .sam file

use List::Util qw(first max maxstr min minstr reduce shuffle sum);

if ($#ARGV ne 2) {
  print "command line: perl bedGraph2sam.pl /seq/chromosome/dm3/dm3.sizes input.bedGraph output.sam\n";
  exit;
}

$tag_size = 25;
$long_fake_seq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";


$chrsize = $ARGV[0];
#$chrdir = $ARGV[1];
$infile = $ARGV[1];
$outfile = $ARGV[2];

$posi_rand = 5;
$strand_rand = 2;

open (in, "<$chrsize");
while ($line=<in>) {
  chomp $line;
  ($chr, $size) = split /\t/, $line;
  $size{$chr} = $size;
  push (@chrs, $chr);
  $num_windows{$chr} = int($size{$chr}/$tag_size);
  print "$chr\t$num_windows{$chr}\n";
}
close in;



$fake_seq = substr($long_fake_seq, 0, $tag_size);

#print "$fake_seq\n";



open (out, ">$outfile");

foreach $chr (@chrs) {

=head
  undef $seq_chr;
  open (in, "<$chrdir/$chr.fa");
  $line = <in>;
  while ($line=<in>) {
    chomp $line;
    $seq_chr = $seq_chr.$line;
  }
  close in;
  print "chromosome $chr loaded\n";
=cut

  undef %flag;
  undef @exp;
    
  my @exp = ();

  open (in, "<$infile");   
  $count = 0;
  while ($line=<in>) {
    chomp $line;
    @data = split /[\s+\t+]/, $line;
    
    if ($data[0] eq $chr) {
      
      $count ++;
      if ($count%10000000 eq 0) {
	print "$chr\t$count\t$line\n";
      }
      
      for ($j=$data[1]; $j<$data[2]; $j++) {
	
	$exp[$j] = $data[3];
	if ($exp[$j] < 0) {
	  $exp[$j] = 0;
	}

      }
    }
  }
  close in;
  
  $count = 0;
       
  for ($i=1; $i<=$num_windows{$chr}; $i++) {
    
    $count ++;
    if ($count%100000 eq 0) {
      print "$chr\t$count\t$line\n";
    }
    
    
    $start = ($i-1)*$tag_size;
    $end = $i*$tag_size;
    
#    print out "$chr\t$start\t$end";
    
    $value = 0;
    @sub_array = @exp[$start..$end];
    $value = sum(@sub_array)/@sub_array;

    if (($start > ($#posi_adj+1)) and ($end < ($size{$chr}-$#posi_adj-1))) {
      for ($j=1; $j<=$value; $j++) {

	$posi_adj = int(rand ($posi_rand));
	$strand_adj = int(rand($strand_rand));
	
	if ($strand_adj eq 1) {
	  $strand = 0;
	  $start_adj = $start + $posi_adj;
	}
	else {
	  $strand = 16;
	  $start_adj = $start - $posi_adj;
	}
	
	$row_counter ++;
	print out "$row_counter\t$strand\t$chr\t$start_adj\t255\t25M\t*\t0\t0\t$fake_seq\t$fake_seq\tXA:i:0\tMD:Z:36 NM:i:0\n";
#	print "$row_counter\t$strand\t$chr\t$start_adj\t255\t25M\t*\t0\t0\t$fake_seq\t$fake_seq\tXA:i:0\tMD:Z:36 NM:i:0\n";

      }
    }
  }

}

close out;
