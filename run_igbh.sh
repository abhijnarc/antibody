#!/bin/bash
# script to run igblast for human database
# define variables
export BLASTDB=/media/AbD/tools/ncbi-igblast-1.22.0/database/human
export IGDATA=/media/AbD/tools/ncbi-igblast-1.22.0/bin
AUXDIR=/media/AbD/tools/ncbi-igblast-1.22.0/optional_file
# get input from command line
inp=$1
# run igblast
${IGDATA}/igblastn -germline_db_V ${BLASTDB}/human_V \
	-germline_db_D ${BLASTDB}/human_D \
	-germline_db_J ${BLASTDB}/human_J \
	-organism human -query $inp \
	-auxiliary_data ${AUXDIR}/human_gl.aux -show_translation -outfmt 3

