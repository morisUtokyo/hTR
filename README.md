# hTR
A pipeline for computing a pair of single nucleotide variants (SNVs) closest to a focal region with tandem repeats of interest. Rather than using publicly available SNV databases, the pipeline attempts to find SNVs among HiFi reads from scratch to identify novel rare SNVs that might be associated with rare TRs.

1. Determine a focal region with tandem repeats in the reference human genome (e.g., hg38). The region can be specified by (chromosome number), (begin), and (end) that respectively represent the human chromosome number, start and end position: e.g., 12,2345,3456 .

2. From reads that are collected from individuals, select reads that map to the focal region, and generate a fasta file with substrings (of length, say 1 kb) of the selected reads before and after the focal region. Each substring of a read must have a annotation of the form "X,readY,[pre|post]" where X is the identifier (alphanumeric string) of an individual, Y is the identifier of the read, and pre (post, respectively) indicates that the substring is located before (after) the region. For example, 
A123,read1,pre, and B456,read4,post.

3. Align each read in the fasta file to the reference genome (e.g., hg38) using minimap2.
- Download the reference human genome hg38 (using say wget_chr.sh), and unzip the files to obtain chromosome DNA sequences named chrZ.fa for each Z=1,2,...,22,X,Y.
- To accelerate computational performance, generate minimap2 index named chrZZZ.mmi on chrZZZ.fa using minimap2_indexing.sh, and put the minimap2 index (chrZZZ.mmi) into a directory, say hg38.
- Generate alignments and put them into a cigar file with cs SAM/PAF tags by calling: 
<br> minimap2 -c --cs chrZZZ.mmi (fasta file)  \>  (cigar file) <br>

4. Execute the Makefile to create an executable named "hap" that parses the cigar file for each substring in the fasta file, and outputs a pair of positions of SNVs that differ from of the reference and are closest to the focal region. Specifically call:　
<br> hap -f (cigarFile) -c (chromosome number) -b (begin) -e (end) \> (hap file) <br>
Each line of (hap file) has the form 
<br> X(tab)Y(tab)H <br> 
where (tab) denotes a tab, X is the ID of an individual, Y is the ID of a read, and H is the pair of positions of closest SNVs separated by a vertical bar "|" : e.g., 
<br> A123   1   2345|3456 .<br> 
When no SNVs are found in a substring before the focal region, H has no position before | : e.g., 
<br> A123  1  |3456 .<br> 
Similarly, in the absence of SNVs before and after the region, H has no positions and is "|".

5. To execute the alignment of substrings in a fasta file to the reference and detection of closest SNVs, call hTR.sh :　
<br> bash hTR.sh (fasta file) (chromosome number) (begin) (end)

