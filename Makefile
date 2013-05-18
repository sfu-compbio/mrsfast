DEBUG := 0
PROFILE := 0
MRSFAST_VERSION := "3.0.0"
BUILD_DATE := "$(shell date)"

all: OPTIMIZE_FLAGS build
debug: DEBUG_FLAGS build
profile: SSE_FLAGS build
build: clean-executable SSE_FLAGS compile mrsfast snp_indexer clean

LDFLAGS=-static
LIBS=-lz -lm -pthread -lpthread
CFLAGS=-DMRSFAST_VERSION=\"$(MRSFAST_VERSION)\" -DBUILD_DATE=\"$(BUILD_DATE)\"

objects=baseFAST.o Sort.o MrsFAST.o Common.o CommandLineParser.o RefGenome.o HashTable.o Reads.o Output.o SNPReader.o HELP.o 

compile: $(objects)

mrsfast:
	gcc -w $(objects) -o $@ ${LDFLAGS} ${LIBS}

snp_indexer: SNPIndexer.o
	gcc $^ -o $@ ${LDFLAGS} ${LIBS}

clean:
	@rm -f $(objects)

clean-executable:
	@rm -f mrsfast

HELP.o:
	@groff -Tascii -man ./HELP.man > HELP
	@ld -r -b binary -o HELP.o HELP

DEBUG_FLAGS:
	$(eval CFLAGS = $(CFLAGS) -ggdb)
	$(eval LIBS = $(LIBS) -ggdb)

OPTIMIZE_FLAGS:
	$(eval CFLAGS = $(CFLAGS) -O2)

PROFILE_FLAGS:
	$(eval CFLAGS = $(CFLAGS) -pg -g)
	$(eval LIBGS = $(LIBS) -pg -g)

SSE_FLAGS:
	$(eval CFLAGS = $(CFLAGS) \
	$(shell gv=`gcc -dumpversion`; \
		sc=`grep -c "sse4" /proc/cpuinfo`; \
		echo $$sc.$$gv | awk -F. '{if($$1>0) print "-DSSE4=1"; else print "-DSSE4=0";  if ($$1>0 && $$2>=4 && $$3>=0) print "-msse4.2" }'))

