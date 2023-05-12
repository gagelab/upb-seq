#!/bin/bash

set -euxo pipefail

# First, make sure you are running this in an environment with the software
#  we need. If this has not been created yet, make it with the .yaml file 
#  Note that this script uses a different environment than the previous few.
#  in the parent directory, e.g.: micromamba create -f freebayes.yaml
# Then, activate with: micromamba activate freebayes


BAMDIR=/mnt/ghowl/upb-seq/align/bams/
REF=/mnt/ghowl/upb-seq/genome/Zm-B73-REFERENCE-NAM-5.0.fa
NPROC=24

# Idex reference genome
if [ ! -f $REF.fai ]
then
  samtools faidx $REF
fi

# Get list of bam filenames
ls -1 $BAMDIR/*.bam > call_snps/bam_list.txt

# Run freebayes
freebayes-parallel \
  <(fasta_generate_regions.py $REF.fai 100000) $NPROC \
  -L call_snps/bam_list.txt \
  -f $REF |\
bgzip \
  -c \
  -@ $NPROC \
> /mnt/ghowl/upb-seq/snps/upb.vcf.gz

