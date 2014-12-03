#!/bin/bash
find ./ -name '*.pdb'  | while read i
do
	echo ${i:2}
	python pdb2fasta.py -t ${i:2}  list.txt
done
