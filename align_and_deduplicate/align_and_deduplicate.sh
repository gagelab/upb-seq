#!/bin/bash

set -euxo pipefail

# In this script, we align each sample's reads
#  against the B73v5 reference genome.
# We will also mark duplicates, ie reads that
#  align to the exact same place and may be
#  PCR duplicates.

# First, make sure you are running this in an environment with the software
#  we need. If this has not been created yet, make it with the .yaml file 
#  in the parent directory, e.g.: micromamba create -f udp_seq.yaml
# Then, activate with: micromamba activate udp-seq

NPROC=24
DATADIR=/mnt/ghowl/upb-seq/
FQDIR=/mnt/ghowl/upb-seq/demultiplex_fastqs/fastqs/
REFPATH=https://download.maizegdb.org/Zm-B73-REFERENCE-NAM-5.0/Zm-B73-REFERENCE-NAM-5.0.fa.gz

# Download and index the B73 reference genome:
REF=$( basename $REFPATH | sed 's/.gz//' )
INDEX=$DATADIR/genome/${REF%%.fa}
if [ ! -f /mnt/ghowl/upb-seq/genome/$REF ]
then
  mkdir -p $DATADIR/genome/
  wget $REFPATH -O - | gzip -dc >  $DATADIR/genome/$REF
  hisat2-build -p $NPROC $DATADIR/genome/$REF $INDEX
fi

# Get list of sample names
SAMPLES=$( find $FQDIR -name *R1.fastq.gz | sed 's/_R1.fastq.gz//' )

mkdir -p $DATADIR/align/bams
mkdir -p $DATADIR/align/logs
for SAMPLE in $SAMPLES
do
  # Don't align unmatched fastqs
  if [ $( basename $SAMPLE ) == "unmatched" ]
  then
    continue
  fi

  # Create name/path for bam file
  OUT=$DATADIR/align/bams/$( basename $SAMPLE ).bam
  # Run hisat2 to align, then sort, fixmates, and mark 
  #  duplicates with samtools
  hisat2 \
    -X 700 \
    --no-spliced-alignment \
    --no-mixed \
    --no-discordant \
    -p $NPROC \
    -x $INDEX \
    -1 ${SAMPLE}_R1.fastq.gz \
    -2 ${SAMPLE}_R2.fastq.gz \
    2>$DATADIR/align/logs/$(basename $SAMPLE).log |\
  awk 'substr($1, 0, 1)=="@" || $0 ~ /NH:i:1/' | \
  samtools sort -@ $NPROC -n |\
  samtools fixmate -@ $NPROC} -m - - |\
  samtools sort -@ $NPROC | \
  samtools markdup -@ $NPROC - - > $OUT
done
