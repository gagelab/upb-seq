#!/bin/bash

DATADIR=/rs1/researchers/j/jlgage/upb-seq
REF=${DATADIR}/genome/Zm-B73-REFERENCE-NAM-5.0.fa
BAMLIST=${DATADIR}/call_snps/bam_list.txt

while read CHUNK
do

  bsub \
    -n 1 \
    -W 48:00 \
    -q gage \
    -o out/$CHUNK.out \
    -e err/$CHUNK.err \
    -J $CHUNK \
  ./call_snps_chunk.sh \
    $REF \
    $BAMLIST \
    $CHUNK \
    ${DATADIR}/vcfs/upb_${CHUNK}.vcf

done < chr_segments.txt
