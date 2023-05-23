#!/bin/bash

set -euxo pipefail

VCFIN=$1
VCFOUT=$2

vcffilter -f "QUAL > 30 & TYPE = snp & NUMALT = 1 & LEN = 1" $VCFIN > $VCFOUT
