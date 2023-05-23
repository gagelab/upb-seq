#!/bin/bash

DATADIR=/rs1/researchers/j/jlgage/upb-seq

for VCF in $DATADIR/vcfs/*.vcf
do

  NAME=$( basename $VCF )
  VCFOUT=$DATADIR/vcfs_filtered/$NAME
  mkdir -p $DATADIR/vcfs_filtered
  bsub \
    -n 1 \
    -W 72:00 \
    -q gage \
    -o out/filter_${NAME%%.vcf}.out \
    -e err/filter_${NAME%%.vcf}.err \
    -J $NAME \
  ./filter_snps_chunk.sh $VCF $VCFOUT

done
