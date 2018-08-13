Hash a fastq or bam in an order independent way.
This allows lossy data transformations to be detected.
Only the name, sequence, and qualities are hashed.
Only the 'primary' reads are hashed.

Raw htslib buffer values are used for hashing, so results may be dependent on version of that library.
The transform to make fastq data equivalent to htslib may therefore become wrong.
Reference to
* https://academic.oup.com/nar/article/38/6/1767/3112533
* http://samtools.github.io/hts-specs/SAMv1.pdf

Problems:
* It does not handle 'U' '.' or '-'
* It does not handle all the non-dna codes

Todo:
* deal with /1 /2 named pairs.
* test with valgrind and debugs
* test O3 vs debug
