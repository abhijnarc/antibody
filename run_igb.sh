#!/bin/bash
# script to run igblast
# define variables
export BLASTDB=/media/AbD/tools/ncbi-igblast-1.22.0/database/mouse
export IGDATA=/media/AbD/tools/ncbi-igblast-1.22.0/bin
AUXDIR=/media/AbD/tools/ncbi-igblast-1.22.0/optional_file
# get input from command line
inp=$1
# run igblast
${IGDATA}/igblastn -germline_db_V ${BLASTDB}/mouse_V \
	-germline_db_D ${BLASTDB}/mouse_D \
	-germline_db_J ${BLASTDB}/mouse_J \
	-organism mouse -query $inp \
	-auxiliary_data ${AUXDIR}/mouse_gl.aux -show_translation -outfmt 19

