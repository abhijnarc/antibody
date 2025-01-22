#!/bin/bash
# script to run igblast in a loop
# define variables
export BLASTDB=/media/AbD/tools/ncbi-igblast-1.22.0/database/mouse
export IGDATA=/media/AbD/tools/ncbi-igblast-1.22.0/bin
AUXDIR=/media/AbD/tools/ncbi-igblast-1.22.0/optional_file
# run loop over director
for infile in `ls fasta/*.fa`
do
	pref=`echo $infile | cut -d/ -f2 | cut -d. -f1`
	outfile="airr/${pref}.airr"
	echo "running igblast for ${pref}...."
	# run igblast
	${IGDATA}/igblastn -germline_db_V ${BLASTDB}/mouse_V \
		-germline_db_D ${BLASTDB}/mouse_D \
		-germline_db_J ${BLASTDB}/mouse_J \
		-organism mouse -query $infile -out $outfile \
		-auxiliary_data ${AUXDIR}/mouse_gl.aux \
		-show_translation -outfmt 19
done
