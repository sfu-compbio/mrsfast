ALL: mrsfast

LDFLAGS=-s -static
LIBS=-lz -lm
CFLAGS= -O2 -s 

mrsfast: baseFAST.o MrsFAST.o Common.o CommandLineParser.o RefGenome.o HashTable.o HashTableBS.o Reads.o Output.o
	gcc $^ -o $@ ${LDFLAGS} ${LIBS}
	rm *.o
