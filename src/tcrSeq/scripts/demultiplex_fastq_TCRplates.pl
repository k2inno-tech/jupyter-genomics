#!/usr/bin/perl
# This script demultiplexes reads multiplexed in a single large fastq file (pair) and saves them into separate fastq files whose names indicate Plate, Well, and TCR isoform (A or B), for instance 01H12B.fastq
# It is called like this: $perl demultiplex_fastq_TCRplates.pl AB-3_S1_L001_R1_001.fastq AB-3_S1_L001_R2_001.fastq
# It will create 2x96 files (one per TCR isoform) per each Plate
# This script will ignore all reads from plates whose code is commented out (see below). This is useful when there is a mixture of TCR genotyping reads and phenotyping reads, for which there is a separate demultiplex script, demultiplex_fastq_phenoplates.pl

$fileR1 = $ARGV[0];
$fileR2 = $ARGV[1];

open(F1,$fileR1);
open(F2,$fileR2);

%plate = (
   "GCAGA" => "01", #uncomment this line if plate code 01 is among the sequences to be demultiplexed
#   "TCGAA" => "02",
#   "AACAA" => "03",
#   "GGTGC" => "04",
#   "TTGGT" => "05",
#   "CATTC" => "06",
#   "ATTGG" => "07",
#   "CGGTT" => "08",
#   "ATCCT" => "09",
#   "ATGTC" => "10",
#   "TCACG" => "11",
#   "AGACC" => "12",
#   "CCCCA" => "13",
#   "GCGCT" => "14",
#   "TCCTT" => "15",
#   "TATAT" => "16",
#   "CGTAA" => "17",
#   "AAGGT" => "18",
#   "AGCTC" => "19",
#   "CTTGC" => "20",
#   "GTATC" => "21",
#   "TATGA" => "22",
#   "CACAC" => "23",
#   "ACACT" => "24",
#   "ACTAC" => "25",
#   "GTTAC" => "26",
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
%TCR = (
   "GTCAC" => "A", # TCRA
   "GAGAT" => "B",
);

foreach $plateID (keys(%plate)) {
   foreach $rowID (keys(%row)) {
      foreach $colID (keys(%col)) {
         foreach $TCRID (keys(%TCR)) {
            $fh = $plate{$plateID}.$row{$rowID}.$col{$colID}.$TCR{$TCRID};
            open $fh, '>', $fh.".fastq"; #open file for writing at the end
         }
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

   $ID = substr($A2, 2, 5); #plate ID barcode
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


   $ID = substr($A2, 9, 5); #row ID barcode
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

   $ID = substr($B2, 7, 5); #TCR ID
#  now find what the true bar code should have been if imperfect match
   $score = 0;
   $trueID = "";
   foreach $key (keys(%TCR)) {
      my $count = ($ID^$key)  =~ tr/\0//;
      if ($count > $score) {      
         $score = $count;
         $trueID = $key
      }
   }
   if ($score >= 4) {#accept $true_plateID as the true plate ID
      $TCRID = $trueID;
   } else {#leave $plateID blank - sequence won't be output
      $TCRID = ""
   }



      
   if (exists $plate{$plateID}  and exists $row{$rowID}  and exists $col{$colID} and exists $TCR{$TCRID}) {
      $fh = $plate{$plateID}.$row{$rowID}.$col{$colID}.$TCR{$TCRID};
      print $fh $A1.$A2.$A3.$A4.$B1.$B2.$B3.$B4;
   };
}
close(F1);
close(F2);
