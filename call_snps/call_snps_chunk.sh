#!/bin/bash

set -euxo pipefail

REF=$1
BAMLIST=$2
CHUNK=$3
OUT=$4

# Make output dir if needed
mkdir -p $( dirname $OUT )

freebayes \
  -L $BAMLIST \
  -f $REF \
  -r $CHUNK |\
bgzip -c |
> $OUT
