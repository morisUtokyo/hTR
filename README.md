# hTR
Computing a pair of SNVs closest to the focal position (with tandem repeats)

Download the reference human genome hg38 (use, say wget_chr.sh) and unzip the files to obtain chrZZZ.fa for each ZZZ=1,2,...,22,X,Y.

Generate minimap2 index on the reference genome hg38 using minimap2_indexing.sh to generate chrZZZ.mmi.

Put the minimap2 index (chrZZZ.mmi) into a directory, say hg38.

Align a set of reads in the input fasta file to the reference using minimap2 to generate alignments in a cigar file with cs SAM/PAF tags

Executing the Makefile creates an executable named "hap" that parses the cigar file and outputs the pair of single nucleotide variants (SNVs) closest to the focal tandem repeat at the input position given in the arguments for each read in the fasta file. Specifically call "hap -f \<cigarFile\> -c \<chromosome number\> -b $\<begin\> -e \<end\>"

To execute the above two steps, use hTR.sh by calling "bash hTR.sh \<seq.fasta\> \<chromosome number\> \<begin\> \<end\>" 

