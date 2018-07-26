#include "xxhash.h"
#include "string.h"
#include "stdio.h"

int main(int argc, char const *argv[])
{
    char* msg1 = "ABCDEFG";
    char* msg2 = "12345";
    XXH64_hash_t a = XXH64(msg1, strlen(msg1), 0);
    XXH64_hash_t b = XXH64(msg2, strlen(msg2), 0);

    printf("%llu %llu\n", a, b);
    return 0;
}
