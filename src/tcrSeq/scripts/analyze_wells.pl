
%plate = (
   "GCAGA" => "01",
#  "TCGAA" => "02",
#  "AACAA" => "03",
#  "GGTGC" => "04",
#  "TTGGT" => "05",
#  "CATTC" => "06",
);
%row = (
   "TAAGC" => "A",
#  "TGCAC" => "B",
#  "CTCAG" => "C",
#  "GGAAT" => "D",
#  "CGAGG" => "E",
#  "AGGAG" => "F",
#  "TGTTG" => "G",
#  "CAACT" => "H",
);
%col = (
#  "GTTCA" => "01",
#  "CAGGA" => "02",
#  "TTATA" => "03",
#  "CCTGT" => "04",
#  "ACCGC" => "05",
   "ACTTA" => "06",
   "GCTAG" => "07",
   "GACGT" => "08",
   "GGCTA" => "09",
#  "GAATG" => "10",
#  "CCAAC" => "11",
#  "GAGAC" => "12",
);
%TCR = (
   "GTCAC" => "A", # TCRA
   "GAGAT" => "B", 
); 


foreach $plateID (sort (keys(%plate))) {
   foreach $rowID (sort (keys(%row))) {
      foreach $colID (sort (keys(%col))) {
         foreach $TCRID (sort (keys(%TCR))) {
            $fh = $plate{$plateID}.$row{$rowID}.$col{$colID}.$TCR{$TCRID};
            print "$fh\n";
            system("java -Xmx10g -jar ./mitcr.jar -pset flex -gene TR$TCR{$TCRID} $fh.fastq $fh\_result.txt")      
         }
      }
   }
}

