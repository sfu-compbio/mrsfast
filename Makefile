ALL: mrsfast

LDFLAGS=-static
LIBS=-lz -lm -g -pg -ggdb
CFLAGS= -O2 -g -pg -ggdb 

mrsfast: baseFAST.o MrsFAST.o Common.o CommandLineParser.o RefGenome.o HashTable.o Reads.o Output.o
	gcc $^ -o $@ ${LDFLAGS} ${LIBS}
	rm *.o

clean:
	rm *.o
