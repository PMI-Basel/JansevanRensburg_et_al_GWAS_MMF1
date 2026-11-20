Sample switch: WTm2 and WTp2
============================

WTp2 and WTm2 were mixed up when pooling the replicate samples for sequencing, so their labels are exchanged here before any downstream analysis.

Operation
---------
swap_samples.R exchanges the column names WTm2 and WTp2 in ReadsPerGene_all_original.txt
and writes ReadsPerGene_all_swapped.txt. Column order is preserved, no values are edited.

Verification
------------
- Only the WTm2 and WTp2 columns differ between the two matrices.
- original$WTm2 is identical to swapped$WTp2, and original$WTp2 to swapped$WTm2.
- Library sizes: WTm2 46.1 M, WTp2 37.9 M assigned reads (unchanged by the swap).
- MD5 of ReadsPerGene_all_swapped.txt: cde4e9685ce3640662dc408fb935f2b8
  - this is the file read by 06_Transcriptome.qmd in the manuscript repository.
- MD5 of ReadsPerGene_all_original.txt: 47c2ef00fd3a5ee98f65323b509c77d2

Consequence
-----------
All results in the manuscript use the swapped (corrected) labelling.
