#!/bin/bsh 
# bash hTR.sh <seq.fasta> <chromosome number> <begin> <end>

fasta_file=$1
chr=$2   # Give a chromosome number, say 1 for chr1
begin=$3 # Give the begin and end pisitions of the focal tandem repeat 
end=$4   # in the 3rd and 4th arguments.

chrFile=../hg38/chr${chr}.mmi
cigarFile=aln.txt
minimap2 -c --cs $chrFile $fasta_file > $cigarFile

hapFile=hap.txt
haplotypeDir=.
${haplotypeDir}/hap -f $cigarFile -c ${chr} -b ${begin} -e ${end} > $hapFile
