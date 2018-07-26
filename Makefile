CPPFLAGS = -g -ggdb3 -I../../download/xxHash -I../../download/htslib
LDLIBS = -lpthread -lz -llzma -lbz2

main: main.o ../../download/xxHash/libxxhash.a ../../download/htslib/libhts.a
