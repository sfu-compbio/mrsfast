GCCMaV := $(shell expr `gcc -dumpversion | cut -f1 -d.` \>= 4)$(shell expr `gcc -dumpversion | cut -f2 -d.` \>= 4)


ALL: mrsfast

LDFLAGS=-static
LIBS=-lz -lm -g -pg -ggdb

CFLAGS= -O2 -g -pg -ggdb 
ifeq "$(GCCMaV)" "11"
	CFLAGS += -msse4.2
endif

mrsfast: baseFAST.o MrsFAST.o Common.o CommandLineParser.o RefGenome.o HashTable.o Reads.o Output.o
	gcc $^ -o $@ ${LDFLAGS} ${LIBS}
	rm *.o

clean:
	echo $(GCCMaV)
	rm *.o
