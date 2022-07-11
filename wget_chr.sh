#!/bin/bash

hg38URL=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes
echo ${hg38URL}

for i in {1..22} ; do
	echo ${i}
	wget ${hg38URL}/chr${i}.fa.gz
done

wget ${hg38URL}/chrX.fa.gz
wget ${hg38URL}/chrY.fa.gz

