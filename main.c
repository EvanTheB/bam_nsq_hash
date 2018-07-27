#include "string.h"
#include "stdio.h"

#include "xxhash.h"
#include "htslib/sam.h"

#ifdef __SANITIZE_ADDRESS__
#include "sanitizer/lsan_interface.h"
#endif

#define bam_get_qname_l(b) ((b)->core.l_qname)
#define bam_get_seq_l(b)  ((b)->core.l_qseq)
#define bam_get_qual_l(b) ((b)->core.l_qseq)
#define bam_is_primary(b) (!((b)->core.flag & (BAM_FSECONDARY|BAM_FSUPPLEMENTARY)))



/*
see seq_comp_table in samtools
import toolz.itertoolz as tz
i = [ 0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15 ]
o = []
for h in range(len(i)):
    for l in range(len(i)):
        o.append((i[l] << 4) + i[h])
print(',\n'.join(', '.join(str(c) for c in l) for l in tz.partition(16, o)))
o = []
for h in range(len(i)):
    for l in range(len(i)):
        o.append((i[h] << 4) + i[l])
print(',\n'.join(', '.join(str(c) for c in l) for l in tz.partition(16, o)))
*/
static uint8_t rcomp1[256] = {
    0, 128, 64, 192, 32, 160, 96, 224, 16, 144, 80, 208, 48, 176, 112, 240,
    8, 136, 72, 200, 40, 168, 104, 232, 24, 152, 88, 216, 56, 184, 120, 248,
    4, 132, 68, 196, 36, 164, 100, 228, 20, 148, 84, 212, 52, 180, 116, 244,
    12, 140, 76, 204, 44, 172, 108, 236, 28, 156, 92, 220, 60, 188, 124, 252,
    2, 130, 66, 194, 34, 162, 98, 226, 18, 146, 82, 210, 50, 178, 114, 242,
    10, 138, 74, 202, 42, 170, 106, 234, 26, 154, 90, 218, 58, 186, 122, 250,
    6, 134, 70, 198, 38, 166, 102, 230, 22, 150, 86, 214, 54, 182, 118, 246,
    14, 142, 78, 206, 46, 174, 110, 238, 30, 158, 94, 222, 62, 190, 126, 254,
    1, 129, 65, 193, 33, 161, 97, 225, 17, 145, 81, 209, 49, 177, 113, 241,
    9, 137, 73, 201, 41, 169, 105, 233, 25, 153, 89, 217, 57, 185, 121, 249,
    5, 133, 69, 197, 37, 165, 101, 229, 21, 149, 85, 213, 53, 181, 117, 245,
    13, 141, 77, 205, 45, 173, 109, 237, 29, 157, 93, 221, 61, 189, 125, 253,
    3, 131, 67, 195, 35, 163, 99, 227, 19, 147, 83, 211, 51, 179, 115, 243,
    11, 139, 75, 203, 43, 171, 107, 235, 27, 155, 91, 219, 59, 187, 123, 251,
    7, 135, 71, 199, 39, 167, 103, 231, 23, 151, 87, 215, 55, 183, 119, 247,
    15, 143, 79, 207, 47, 175, 111, 239, 31, 159, 95, 223, 63, 191, 127, 255
};

static uint8_t comp1[256] = {
    0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15,
    128, 136, 132, 140, 130, 138, 134, 142, 129, 137, 133, 141, 131, 139, 135, 143,
    64, 72, 68, 76, 66, 74, 70, 78, 65, 73, 69, 77, 67, 75, 71, 79,
    192, 200, 196, 204, 194, 202, 198, 206, 193, 201, 197, 205, 195, 203, 199, 207,
    32, 40, 36, 44, 34, 42, 38, 46, 33, 41, 37, 45, 35, 43, 39, 47,
    160, 168, 164, 172, 162, 170, 166, 174, 161, 169, 165, 173, 163, 171, 167, 175,
    96, 104, 100, 108, 98, 106, 102, 110, 97, 105, 101, 109, 99, 107, 103, 111,
    224, 232, 228, 236, 226, 234, 230, 238, 225, 233, 229, 237, 227, 235, 231, 239,
    16, 24, 20, 28, 18, 26, 22, 30, 17, 25, 21, 29, 19, 27, 23, 31,
    144, 152, 148, 156, 146, 154, 150, 158, 145, 153, 149, 157, 147, 155, 151, 159,
    80, 88, 84, 92, 82, 90, 86, 94, 81, 89, 85, 93, 83, 91, 87, 95,
    208, 216, 212, 220, 210, 218, 214, 222, 209, 217, 213, 221, 211, 219, 215, 223,
    48, 56, 52, 60, 50, 58, 54, 62, 49, 57, 53, 61, 51, 59, 55, 63,
    176, 184, 180, 188, 178, 186, 182, 190, 177, 185, 181, 189, 179, 187, 183, 191,
    112, 120, 116, 124, 114, 122, 118, 126, 113, 121, 117, 125, 115, 123, 119, 127,
    240, 248, 244, 252, 242, 250, 246, 254, 241, 249, 245, 253, 243, 251, 247, 255
};

static void reverse(uint8_t *str, int32_t len)
{
    int32_t high = len - 1;
    int32_t low = 0;
    uint8_t t;
    for (; high > low; low++, high--)
    {
        t = str[high];
        str[high] = str[low];
        str[low] = t;
    }
}

static void nibble_fiddle(uint8_t *seq, int32_t len)
{
    for (int32_t i = 1; i < len; ++i)
    {
        seq[i-1] = (seq[i-1] & 0xF0) | (seq[i] & 0xF);
    }
    seq[len-1] &= 0xF0;
}

static void reverse_complement(uint8_t *seq, int32_t len)
{
    int32_t byte_len = (len + 1) / 2;
    /* reverse seq  */
    reverse(seq, byte_len);

    /* account for odd length  */
    if (len & 1)
    {
        nibble_fiddle(seq, byte_len);
    }

    /* reverse? complement the nibbles
       no reverse on odd length sequence, fiddle did it.
    */
    if (len & 1)
    {
        for (int i = 0; i < byte_len; ++i)
        {
            seq[i] = comp1[seq[i]];
        }
    } else
    {
        for (int i = 0; i < byte_len; ++i)
        {
            seq[i] = rcomp1[seq[i]];
        }
    }
}

int main(int argc, char const *argv[])
{

    samFile *f = sam_open(argv[1], "r");

    bam_hdr_t *h = sam_hdr_read(f);

    bam1_t *b1 = bam_init1();

    XXH64_hash_t acc = 0;
    while (sam_read1(f, h, b1) >= 0)
    {
        if (!bam_is_primary(b1))
        {
            continue;
        }

        if (bam_is_rev(b1))
        {
            reverse_complement(bam_get_seq(b1), bam_get_seq_l(b1));
            reverse(bam_get_qual(b1), bam_get_qual_l(b1));
        }

        acc += XXH64(bam_get_qname(b1), bam_get_qname_l(b1), 0);
        acc += XXH64(bam_get_seq(b1), (bam_get_seq_l(b1) + 1) / 2, 0);
        acc += XXH64(bam_get_qual(b1), bam_get_qual_l(b1), 0);
    }

    bam_destroy1(b1);
    bam_hdr_destroy(h);
    sam_close(f);

    printf("%llx\n", acc);

    #ifdef __SANITIZE_ADDRESS__
    __lsan_do_leak_check();
    __lsan_disable();
    #endif

    return 0;
}
