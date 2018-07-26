#include "string.h"
#include "stdio.h"

#include "xxhash.h"
#include "htslib/sam.h"

#define bam_get_qname_l(b) ((b)->core.l_qname)
#define bam_get_seq_l(b) (((b)->core.l_qseq + 1) / 2)
#define bam_get_qual_l(b) ((b)->core.l_qseq)

int main(int argc, char const *argv[])
{

    samFile *f = sam_open(argv[1], "r");

    bam_hdr_t *h = sam_hdr_read(f);

    bam1_t *b1 = bam_init1();

    XXH64_hash_t acc = 0;
    while (sam_read1(f, h, b1) >= 0)
    {
        acc += XXH64(bam_get_qname(b1), bam_get_qname_l(b1), 0);
        acc += XXH64(bam_get_seq(b1), bam_get_seq_l(b1), 0);
        acc += XXH64(bam_get_qual(b1), bam_get_qual_l(b1), 0);
    }

    printf("%llx\n", acc);
    return 0;
}
