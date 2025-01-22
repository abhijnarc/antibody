#!/bin/bash

# try to generate file names
for infile in `ls fasta/heavy_*.fa`
do
	contig=`echo $infile | cut -d/ -f2 | cut -d. -f1`
	echo $contig
done
