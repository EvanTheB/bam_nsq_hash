CFLAGS = -g -Wall -std=gnu99
LDLIBS = -lhts -lxxhash -lpthread -lz -llzma -lbz2 -ldeflate -lcurl -lcrypto

release: CFLAGS += -O3
release: bam_nsq_hash

bam_nsq_hash:

debug: CFLAGS += -ggdb3 -fsanitize=address -fsanitize=undefined
debug: LDFLAGS += -fsanitize=address -fsanitize=undefined
debug: bam_nsq_hash

pedantic: CFLAGS += -Wextra -pedantic -std=c99
# pedantic: CC = clang
# pedantic: CFLAGS += -Wextra -Weverything -pedantic -std=c99
pedantic: bam_nsq_hash

.PHONY: clean test
clean:
	rm bam_nsq_hash

test: bam_nsq_hash
	cd test && bash test.sh ../bam_nsq_hash
