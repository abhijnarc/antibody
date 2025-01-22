#!/bin/bash

# fastp for a single sample
# get prefix from command line
samp=$1
echo "Sample .... $samp"
r1_file="${samp}_R1_001.fastq.gz"
r2_file="${samp}_R2_001.fastq.gz"
o1_file="${samp}_fastp_R1_001.fastq.gz"
o2_file="${samp}_fastp_R2_001.fastq.gz"
rep_file="${samp}_fastp.html"
/home/user/tools/fastp --in1 $r1_file --in2 $r2_file \
		--out1 $o1_file --out2 $o2_file \
		--qualified_quality_phred 30 \
		--detect_adapter_for_pe \
		--correction \
		--trim_tail1=1 \
		--cut_tail \
		--cut_window_size=4 \
		--cut_mean_quality=30 \
		--length_required=50 \
		--html $rep_file
