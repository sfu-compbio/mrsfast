GCC44 := $(shell expr `gcc -dumpversion | cut -f1 -d.` \>= 4)$(shell expr `gcc -dumpversion | cut -f2 -d.` \>= 4)
SSE4 := $(shell expr "`cat /proc/cpuinfo | grep sse4`" != "")
DEBUG := 0
PROFILE := 0
MRSFAST_VERSION := "3.0.0"
BUILD_DATE := "$(shell date)"

ALL: mrsfast snp_indexer clean

LDFLAGS=-static
LIBS=-lz -lm -pthread -lpthread
CFLAGS=-O2 -DSSE4=$(SSE4) -DMRSFAST_VERSION=\"$(MRSFAST_VERSION)\" -DBUILD_DATE=\"$(BUILD_DATE)\"


ifeq "$(GCC44)$(SSE4)" "111"
	CFLAGS += -msse4.2
endif

ifeq "$(DEBUG)" "1"
	CFLAGS += -ggdb
	LIBS += -ggdb
endif

ifeq "$(PROFILE)" "1"
	CFLAGS += -pg
	LIBS += -pg
endif

mrsfast: baseFAST.o Sort.o MrsFAST.o Common.o CommandLineParser.o RefGenome.o HashTable.o Reads.o Output.o SNPReader.o HELP.o 
	gcc $^ -o $@ ${LDFLAGS} ${LIBS}

snp_indexer: SNPIndexer.o
	gcc $^ -o $@ ${LDFLAGS} ${LIBS}

clean:
	rm -rf *.o

HELP.o:
	ld -r -b binary -o HELP.o HELP
