MRSFAST_VERSION := "3.4.6"
BUILD_DATE := "$(shell date)"

all: OPTIMIZE_FLAGS build
debug: DEBUG_FLAGS build
profile: PROFILE_FLAGS build
build: clean_executables mrsfast snp_indexer clean_objects

LDFLAGS=#-static
LIBS=-lz -lm -pthread -lpthread  -DSSE4=1 -msse4.2
CFLAGS=-DMRSFAST_VERSION=\"$(MRSFAST_VERSION)\" -DBUILD_DATE=\"$(BUILD_DATE)\" -DSSE4=1 -msse4.2

objects=baseFAST.o Sort.o MrsFAST.o Common.o CommandLineParser.o RefGenome.o HashTable.o Reads.o Output.o SNPReader.o  

mrsfast: clean_executables $(objects)
ifeq ($(shell uname -s),Linux)
	gcc -w $(objects) -o $@ ${LDFLAGS} ${LIBS}
else
	gcc -Wl,-no_pie -fno-pic -w $(objects) -o $@ ${LDFLAGS} ${LIBS}
endif

snp_indexer: clean_executables SNPIndexer.o
	gcc SNPIndexer.o -o $@ ${LDFLAGS} ${LIBS}

clean_objects: mrsfast snp_indexer
	@rm -f $(objects)
	@rm -f SNPIndexer.o
	@rm -f HELPstub.c
	@rm -f HELPstub.o

clean:
	@rm -f $(objects)
	@rm -f SNPIndexer.o
	@rm -f HELPstub.c
	@rm -f HELPstub.o
	
clean_executables:
	@rm -f mrsfast
	@rm -f snp_indexer


DEBUG_FLAGS:
	$(eval CFLAGS = $(CFLAGS) -ggdb)
	$(eval LIBS = $(LIBS) -ggdb)

OPTIMIZE_FLAGS:
	$(eval CFLAGS = $(CFLAGS) -O2)

PROFILE_FLAGS:
		$(eval CFLAGS = $(CFLAGS) -pg -g)
	$(eval LIBS = $(LIBS) -pg -g)

