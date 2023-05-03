#!/bin/bash

set -euxo pipefail

# Generally following guidelines from Twist outlined here:
#  https://www.twistbioscience.com/sites/default/files/resources/2022-01/Guide_NGS_96Plex_LibraryPrepKit_SampleDemultiplexingGiode_24JAN22_Rev1.0.pdf

# These sequencing libraries are multiplexed in a manner that does not get
#  demultiplexed by the standard procedures at GSL, so we need to do it 
#  manually.

# First, make sure you are running this in an environment with the software
#  we need. If this has not been created yet, make it with the .yaml file 
#  in this directory, e.g.: micromamba create -f demultiplex_fastqs.yaml
# Then, activate with: micromamba activate demultiplex_fastqs

# Can set these variables if running on a different computer or want
#  to output results to a different directory.
RAWFQDIR=/mnt/ghowl/RawSeqData/zea/wgs/UPB_seq/NVS174A_L[34]_RellenAlvarez/
SAMPLESHEETDIR=./sample_sheets/
SAMPLESHEETBASE=tsi_seq_sample_sheet_
OUTDIR=/mnt/ghowl/upb-seq/demultiplex_fastqs/
MAXMISMATCH=1

# This will only find .fastq.gz files, and finds the unique library names
#  without R1, R2, etc.
LIBNAMES=$( ls -1 $RAWFQDIR/*R1*fastq.gz | sed 's/_R1.*fastq.gz//' | sort | uniq )

# Make directories for output
mkdir -p $OUTDIR/fastqs
mkdir -p $OUTDIR/metrics

# Loop through libraries
for LIB in $LIBNAMES
do

  # Get paths to read1 and read2
  R1=$( ls ${LIB}*R1* )
  R2=$( ls ${LIB}*R2* )

  # Pull the library number out from the library name
  SAMPLENUM=$( echo $LIB | grep -Po "Library\K[0-9]" )
  # Get the path to the sample sheet
  SAMPLESHEET=$SAMPLESHEETDIR/${SAMPLESHEETBASE}$SAMPLENUM.csv

  # Set paths for where to put demultiplexed fastqs...
  OUTPREFIX=${OUTDIR}/fastqs/$( basename $LIB )
  # ... and metrix about the dumultiplexing
  OUTMETRICS=${OUTDIR}/metrics/$( basename $LIB )

  # Run demultiplexing software
  fgbio -Xmx32G DemuxFastqs \
    --inputs $R1 $R2 \
    --metadata $SAMPLESHEET \
    --read-structures 8B12S+T 8S+T \
    --max-mismatches=$MAXMISMATCH \
    --output $OUTPREFIX \
    --metrics $OUTMETRICS \
    --threads $(nproc) \
    --output-type Fastq
  
done
