CFLAGS = -g -Wall -std=gnu99
CPPFLAGS = -I../../download/xxHash -I../../download/htslib
LDLIBS = -lpthread -lz -llzma -lbz2

main: main.o ../../download/xxHash/libxxhash.a ../../download/htslib/libhts.a

release: CFLAGS += -O2
release: main

debug: CFLAGS += -ggdb3 -fsanitize=address -fsanitize=undefined
debug: LDFLAGS += -fsanitize=address -fsanitize=undefined
debug: main

clean:
	rm main.o main
