CFLAGS ?= -g -ggdb3
CPPFLAGS ?= -I../../download/xxHash -I../../download/htslib
LDLIBS ?= -lpthread -lz -llzma -lbz2

main: main.o ../../download/xxHash/libxxhash.a ../../download/htslib/libhts.a
