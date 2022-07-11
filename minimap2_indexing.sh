#!/bin/bash

for i in {1..22} ; do
	echo ${i}
	minimap2 -d chr${i}.mmi chr${i}.fa
done

minimap2 -d chrX.mmi chrX.fa
minimap2 -d chrY.mmi chrY.fa
