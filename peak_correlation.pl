#!/usr/bin/perl

=head

this script is used to find chip-seq correlation within a peak range

input files are: chip-seq_BedGraph_file

=cut

use List::Util qw(first max maxstr min minstr reduce shuffle sum);


if ($#ARGV ne 3) {
  print "command line: perl peak_correlation.pl peak.xls even.bedGraph odd.bedGrap merge.bedGraph\n";
  exit;
}


$peakfile = $ARGV[0];
$chip_file1 = $ARGV[1];
$chip_file2 = $ARGV[2];
$chip_file3 = $ARGV[3];
$outfile = "$peakfile.corr.xls";


$window_width_thrd = 2000;
$res = 50;
#$num_windows = $half_width/$res;
$flag_log = 1;

# dividing wiggle files to different chromosome
open (in, "<$peakfile");
while ($line=<in>) {
  if ($line=~/^\#/) {
  }
  elsif ($line eq "\n") {
  }
  elsif ($line=~/start/) {
    chomp $line;
    $header = $line;
  }
  else {
    chomp $line;
    @data = split /[\t+\s+]/, $line;
    $chr{$data[0]} = 1;
  }
}
close in; 
@chrs = keys %chr;

print "@chrs\n";

open (out, ">$outfile");
print out "$header\tcorrelation\tfold_differences\taver_coverage\n";

#@chrs = ("chr2L");

foreach $chr (@chrs) {

  undef %flag;
  undef @exp1;
  undef @exp2;
  undef @exp3;

  my @exp1 = ();
  my @exp2 = ();
  my @exp3 = ();

  print "$chr\n";

  open (in, "<$peakfile");
  while ($line=<in>) {
    if (!(($line=~/^\#/) or ($line=~/start/))) {
      chomp $line;
      @data = split /[\t+\s+]/, $line;
      if ($data[0] eq $chr) {
	$start = $data[1] - $window_width_thrd;
	$end = $data[2] + $window_width_thrd;
	for ($i=$start; $i<=$end; $i++) {
	  $flag{$i} = 1;
	}
      }
    }
  }
  close in; 


  open (in, "<$chip_file1");   
  $count = 0;
  while ($line=<in>) {
    chomp $line;
    @data = split /[\s+\t+]/, $line;
    
    if ($data[0] eq $chr) {
      
      $count ++;
      if ($count%100000 eq 0) {
	print "$chr\t$count\t$line\n";
      }
      
      for ($j=$data[1]; $j<$data[2]; $j++) {
	
	if ($flag{$j}) {
	  $exp1[$j] = $data[3];
	  if ($exp1[$j] < 0) {
	    $exp1[$j] = 0;
	  }
	}
      }      
    }
  }
  close in;
  
  open (in, "<$chip_file2");   
  $count = 0;
  while ($line=<in>) {
    chomp $line;
    @data = split /[\s+\t+]/, $line;
    
    if ($data[0] eq $chr) {
      
      $count ++;
      if ($count%100000 eq 0) {
	print "$chr\t$count\t$line\n";
      }
      
      for ($j=$data[1]; $j<$data[2]; $j++) {
	
	if ($flag{$j}) {
	  $exp2[$j] = $data[3];
	  if ($exp2[$j] < 0) {
	    $exp2[$j] = 0;
	  }
	}
      }      
    }
  }
  close in;
  
  open (in, "<$chip_file3");   
  $count = 0;
  while ($line=<in>) {
    chomp $line;
    @data = split /[\s+\t+]/, $line;
    
    if ($data[0] eq $chr) {
      
      $count ++;
      if ($count%100000 eq 0) {
	print "$chr\t$count\t$line\n";
      }
      
      for ($j=$data[1]; $j<$data[2]; $j++) {
	
	if ($flag{$j}) {
	  $exp3[$j] = $data[3];
	  if ($exp3[$j] < 0) {
	    $exp3[$j] = 0;
	  }
	}
      }      
    }
  }
  close in;  


  open (in, "<$peakfile");
  while ($line=<in>) {
    if (!(($line=~/^\#/) or ($line=~/start/))) {
      
      @x = ();
      @y = ();

      chomp $line;
      @data = split /[\t+\s+]/, $line;

      if ($data[0] eq $chr) {
	
	$row_count ++;
	if ($row_count%100 eq 0) {
	  print "$row_count\n";
	}

	$peak_start = $data[1];
	$peak_end = $data[2];
	$peak_length = $data[3];
	$peak_summit = $data[1] + $data[4];
	
	if ($peak_length > $window_width_thrd) {
	  $start = $peak_summit - $window_width_thrd;
	  $num_windows = int(2*$window_width_thrd/$res);
	  $window_size = 2*$window_width_thrd + 1;
	}
	else {
	  $start = $peak_start;
	  $num_windows = int($peak_length/$res);
	  $window_size = $peak_length + 1;
	}

	$s_x = 0;
	$s_y = 0;
	$s_z = 0;

	for ($i=1; $i<=$num_windows; $i++) {
    
	  $w_start = ($i-1)*$res + 1 + $start;
	  $w_end = $i*$res + $start;
     
	  $value1 = 0;
	  $value2 = 0;
	  $value3 = 0;	 

	  @sub_array1 = @exp1[$w_start..$w_end];
	  @sub_array2 = @exp2[$w_start..$w_end];
	  @sub_array3 = @exp3[$w_start..$w_end];

	  $value1 = sum(@sub_array1);
	  $value2 = sum(@sub_array2);
	  $value3 = sum(@sub_array3);

	  if ($value1<=1) {
	    $value1 = 1;
	  }
	  if ($value2<=1) {
	    $value2 = 1;
	  }
	  if ($value3<=1) {
	    $value3 = 1;
	  }

	  $s_x = $s_x + $value1;
	  $s_y = $s_y + $value2;
	  $s_z = $s_z + $value3;
	 
	  if ($flag_log) {
	    $value1 = log($value1)/log(2);
	    $value2 = log($value2)/log(2);
	  }    
	  if ($value1 or $value2) {
	    push (@x, $value1);
	    push (@y, $value2);
	  }

	}

	$aver_coverage = $s_z/$window_size;
#	$s_x = sum(@x);
#	$s_y = sum(@y);

	if ($s_x and $s_y) {
	  if ($s_x>=$s_y) {
	    $s_ratio = $s_x/$s_y;
	  }
	  else {
	    $s_ratio = $s_y/$s_x;
	  }
	}
	else {
	  $s_ratio = "inf";
	}

#	print "$x\n$y\n";

	&pearson_correlation;
	print out "$line\t$corr\t$s_ratio\t$aver_coverage\n";
      } 
    }
  }
}

close out;



sub pearson_correlation {
  
  if ($#x ne $#y) {
    print "array x and y are not equal size\n";
    last;
  }
  
  my $n = $#x + 1;
  my $sum_x = 0;
  my $sum_y = 0;
  my $sum_xy = 0;
  my $sum_yy = 0;
  my $sum_xx = 0;

  for (my $i=0; $i<=$#x; $i++) {
    $sum_x = $sum_x + $x[$i];
    $sum_y = $sum_y + $y[$i];
    $sum_xy = $sum_xy + $x[$i]*$y[$i];
    $sum_xx = $sum_xx + $x[$i]*$x[$i];
    $sum_yy = $sum_yy + $y[$i]*$y[$i];
  }

  if ((($n*$sum_xx-$sum_x*$sum_x)>0) and (($n*$sum_yy-$sum_y*$sum_y)>0)) {
    $corr = ($n*$sum_xy-$sum_x*$sum_y)/sqrt($n*$sum_xx-$sum_x*$sum_x)/sqrt($n*$sum_yy-$sum_y*$sum_y);
  }
  else {
    $corr = 0;
  }
  return;

}
