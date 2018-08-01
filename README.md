Hash a fastq or bam in an order independent way.
This allows lossy data transformations to be detected.
Only the name, sequence, and qualities are hashed.
Only the 'primary' reads are hashed.

Raw htslib buffer values are used for hashing, so values may be dependent on version.
The transform to make fastq data equivalent to htslib may also be dependent on version.

known issue: odd length reads have an extra nibble (0), this is the same as a '='.
    This should be fine as the quality has different length anyway.
    Another workaround would be to hash the length as part of the data.

Todo:
    deal with /1 /2 named pairs.
    probably should hash the actual sequence not the htslib nibbles
