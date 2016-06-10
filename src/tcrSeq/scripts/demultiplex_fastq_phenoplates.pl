#!/usr/bin/perl
# This script demultiplexes reads multiplexed in a single large fastq file (pair) and saves them into separate fastq files whose names indicate Plate, Well, and R1 or R2, for instance 01H12R1.fastq
# It is called like this: $perl demultiplex_fastq_phenoplates.pl AB-3_S1_L001_R1_001.fastq AB-3_S1_L001_R2_001.fastq
# It will create 2x96 files per each Plate
# This script will ignore all reads from plates whose code is commented out (see below). This is useful when there is a mixture of TCR genotyping reads and phenotyping reads.

$fileR1 = $ARGV[0];
$fileR2 = $ARGV[1];

open(F1,$fileR1);
open(F2,$fileR2);

%plate = (
#   "GCAGA" => "01",
#   "TCGAA" => "02",
   "AACAA" => "03",
#   "GGTGC" => "04",
#   "TTGGT" => "05",
#   "CATTC" => "06",
);
%row = (
   "TAAGC" => "A",
   "TGCAC" => "B",
   "CTCAG" => "C",
   "GGAAT" => "D",
   "CGAGG" => "E",
   "AGGAG" => "F",
   "TGTTG" => "G",
   "CAACT" => "H",
);
%col = (
   "GTTCA" => "01",
   "CAGGA" => "02",
   "TTATA" => "03",
   "CCTGT" => "04",
   "ACCGC" => "05",
   "ACTTA" => "06",
   "GCTAG" => "07",
   "GACGT" => "08",
   "GGCTA" => "09",
   "GAATG" => "10",
   "CCAAC" => "11",
   "GAGAC" => "12",
);

foreach $plateID (keys(%plate)) {
   foreach $rowID (keys(%row)) {
      foreach $colID (keys(%col)) {
            $fh = $plate{$plateID}.$row{$rowID}.$col{$colID};
            $fh1 = $plate{$plateID}.$row{$rowID}.$col{$colID}."1";
            $fh2 = $plate{$plateID}.$row{$rowID}.$col{$colID}."2";
            open $fh1, '>', $fh."R1.fastq"; #open file for writing at the end
            open $fh2, '>', $fh."R2.fastq"; #open file for writing at the end
      }
   }
}

while($A1 = <F1>) { #read 4 lines from R1 and 4 lines from R2
   $A2 = <F1>;
   $A3 = <F1>;
   $A4 = <F1>;

   $B1 = <F2>;
   $B2 = <F2>;
   $B3 = <F2>;
   $B4 = <F2>;

#  now find out if the bar codes make sense

   $ID = substr($A2, 2, 5); #plate ID
#  now find what the true bar code should have been if imperfect match
   $score = 0;
   $trueID = "";
   foreach $key (keys(%plate)) {
      my $count = ($ID^$key)  =~ tr/\0//;
      if ($count > $score) { 
         $score = $count;
         $trueID = $key
      }
   }
   if ($score >= 4) {#accept $true_plateID as the true plate ID
      $plateID = $trueID;
   } else {#leave $plateID blank - sequence won't be output
      $plateID = ""
   }


   $ID = substr($A2, 9, 5); #row ID
#  now find what the true bar code should have been if imperfect match
   $score = 0;
   $trueID = "";
   foreach $key (keys(%row)) {
      my $count = ($ID^$key)  =~ tr/\0//;
      if ($count > $score) {      
         $score = $count;
         $trueID = $key
      }
   }
   if ($score >= 4) {#accept $true_plateID as the true plate ID
      $rowID = $trueID;
   } else {#leave $plateID blank - sequence won't be output
      $rowID = ""
   }

   $ID = substr($B2, 2, 5); #column ID
#  now find what the true bar code should have been if imperfect match
   $score = 0;
   $trueID = "";
   foreach $key (keys(%col)) {
      my $count = ($ID^$key)  =~ tr/\0//;
      if ($count > $score) {      
         $score = $count;
         $trueID = $key
      }
   }
   if ($score >= 4) {#accept $true_plateID as the true plate ID
      $colID = $trueID;
   } else {#leave $plateID blank - sequence won't be output
      $colID = ""
   }

   if (exists $plate{$plateID}  and exists $row{$rowID}  and exists $col{$colID} ) {
      $fh1 = $plate{$plateID}.$row{$rowID}.$col{$colID}."1";
      $fh2 = $plate{$plateID}.$row{$rowID}.$col{$colID}."2";
      print $fh1 $A1.$A2.$A3.$A4;
      print $fh2 $B1.$B2.$B3.$B4;
   };
}
close(F1);
close(F2);
