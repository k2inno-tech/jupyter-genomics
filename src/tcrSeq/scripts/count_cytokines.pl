%cyt = (
   "GCCGGAGGAGGTGGATGTGC" => "GATA3",
   "CCCAACACAGGAGCGCACTG" => "TBET",
   "GGCAGCCAAGGCCCTGTCGT" => "FOXP3",
   "AGAGGAAGTCCATGTGGGAG" => "RORC",
   "GCGAGCTGGTGCGCACCGAC" => "RUNX1",
   "GGACCACGCAGGCGAGCTCG" => "RUNX3",
   "CCTACACGGCCCCACCTGCC" => "BCL6",
   "CCACAGAACTGAAACATCTT" => "IL2",
   "CCCAAGCTGAGAACCAAGAC" => "IL10",
   "AGACCTCTTTTATGATGGCC" => "IL12A",
   "GGTATGGAGCATCAACCTGA" => "IL13",
   "CAACCTGAACATCCATAACC" => "IL17A",
   "GGGTTCTCTTGGCTGTTACT" => "IFNG",
   "GGAGGCGCTCCCCAAGAAGA" => "TNFA",
   "CCGAGAAGCGGTACCTGAAC" => "TGFB",
   "GCCAACTTTGCAGCCCAGAA" => "PRF1",
   "CCACAATATCAAAGAACAGG" => "GZMB",
);

#opendir(THISDIR,".");
#@allfiles = grep -T, readdir THISDIR;
#$n = scalar(@allfiles);

%plate = (
#  "GCAGA" => "01",
#  "TCGAA" => "02",
   "AACAA" => "03",
#  "GGTGC" => "04",
#  "TTGGT" => "05",
#  "CATTC" => "06",
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

foreach $plateID (sort (keys(%plate))) {
   foreach $rowID (sort (keys(%row))) {
      foreach $colID (sort (keys(%col))) {
         $fh = $plate{$plateID}.$row{$rowID}.$col{$colID};
         open(F1,$fh."R1.fastq");
         open $fh, '>', $fh."R1.count"; #open file for writing at the end
         print $fh "\t$fh\n"; #print header

#        zero out counters
         foreach $key (keys(%cyt)) {$count{$cyt{$key}} = 0};
         while($A1 = <F1>) { #read 4 lines from R1 and 4 lines from R2
            $A2 = <F1>;
            $A3 = <F1>;
            $A4 = <F1>;
#           now find out if the bar codes make sense
            $seq = substr($A2, 36, 20);
            if (exists $cyt{$seq}) {$count{$cyt{$seq}}++}; #add to count
         };
         foreach $key (keys(%cyt)) {
            print $fh $cyt{$key}."\t".$count{$cyt{$key}}."\n"
         };
         close(F1);
         close($fh);

      }
   }
}

#   }
#}
