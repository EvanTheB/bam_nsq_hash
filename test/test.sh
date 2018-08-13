#!/bin/bash
set -euo pipefail
# set -x

DUT=$1
echo "testing" $DUT

# test that it works (+ve)
REF=$(cat na12878.sam.hash)
test $($DUT na12878.sam) = $REF
samtools fastq -n na12878.sam > na12878.fq
test $($DUT -f na12878.fq) = $REF
test $(gzip -c na12878.fq | $DUT -f -) = $REF
test $(bgzip -c na12878.fq | $DUT -f -) = $REF
test $(samtools view na12878.sam -b | $DUT -) = $REF
# test $(samtools view na12878.sam -C | $DUT -) = $REF
test $(samtools sort -n na12878.sam | $DUT -) = $REF

# test that it works (-ve)
test $(head -n-1 na12878.sam | $DUT -) != $REF
test $(sed 's|35203|35204|g' na12878.sam | $DUT -) != $REF
test $(sed 's|AGAGTG|AGAGTC|g' na12878.sam | $DUT -) != $REF
test $(sed 's|0@:=|1@:=|g' na12878.sam | $DUT -) != $REF

# test for difference between different length samples
# even if the longer has value '0' (=) for the last nibble
test $($DUT odd_eq.sam) != $($DUT odd_eq_2.sam)

# check that reverse complement works
test $($DUT og.sam) = $($DUT rcompl.sam)

# a hard fastq
test $($DUT -f hard.fastq) = $(fastqtobam qualitymax=94 hard.fastq | $DUT -)

# upper lower
REF=$($DUT -f upper.fastq)
test $(fastqtobam upper.fastq | $DUT -) = $REF
test $($DUT -f lower.fastq) = $REF
test $(fastqtobam lower.fastq | $DUT -) = $REF
