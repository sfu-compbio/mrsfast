DEBUG := 0
PROFILE := 0
MRSFAST_VERSION := "3.3.6"
BUILD_DATE := "$(shell date)"

all: OPTIMIZE_FLAGS build
debug: DEBUG_FLAGS build
profile: PROFILE_FLAGS build
build: clean-executable SSE_FLAGS compile mrsfast snp_indexer clean

LDFLAGS=#-static
LIBS=-lz -lm -pthread -lpthread
CFLAGS=-fno-pic -DMRSFAST_VERSION=\"$(MRSFAST_VERSION)\" -DBUILD_DATE=\"$(BUILD_DATE)\"

objects=baseFAST.o Sort.o MrsFAST.o Common.o CommandLineParser.o RefGenome.o HashTable.o Reads.o Output.o SNPReader.o  HELP.o

compile: $(objects)

mrsfast:
ifeq ($(shell uname -s),Linux)
	gcc -w $(objects) -o $@ ${LDFLAGS} ${LIBS}
else
	gcc -Wl,-no_pie -fno-pic -w $(objects) -o $@ ${LDFLAGS} ${LIBS}
endif

snp_indexer: SNPIndexer.o
	gcc $^ -o $@ ${LDFLAGS} ${LIBS}

clean:
	@rm -f $(objects)
	@rm -f SNPIndexer.o
	@rm -f HELPstub.c
	@rm -f HELPstub.o
	
clean-executable:
	@rm -f mrsfast

HELP.o:
	@groff -Tascii -man ./HELP.man > HELP
ifeq ($(shell uname -s),Linux)
	@ld -r -b binary -o HELP.o HELP
else
	@touch HELPstub.c
	gcc -o HELPstub.o -c HELPstub.c
	ld -r -o HELP.o -sectcreate binary HELP HELP HELPstub.o
endif

DEBUG_FLAGS:
	$(eval CFLAGS = $(CFLAGS) -ggdb)
	$(eval LIBS = $(LIBS) -ggdb)

OPTIMIZE_FLAGS:
	$(eval CFLAGS = $(CFLAGS) -O2)

PROFILE_FLAGS:
		$(eval CFLAGS = $(CFLAGS) -pg -g)
	$(eval LIBS = $(LIBS) -pg -g)

SSE_FLAGS:
ifeq ($(shell uname -s),Linux)
ifeq ($(with-sse4),no)
		$(shell echo "-DSSE4=0")
else
        	$(eval CFLAGS = $(CFLAGS) \
        	$(shell gv=`gcc -dumpversion`; \
            	    sc=`grep -c "sse4" /proc/cpuinfo`; \
                	echo $$sc.$$gv | awk -F. '{if($$1>0 && $$2>=4 && $$3>=4) print "-DSSE4=1 -msse4.2"; else print "-DSSE4=0"}'))
endif
else
ifeq ($(with-sse4),no)
		$(shell echo "-DSSE4=0")
else
        $(eval CFLAGS = $(CFLAGS) \
        $(shell gv=`gcc -dumpversion`; \
                sc=`sysctl -n machdep.cpu.features | grep -c "SSE4"` ;\
                echo $$sc.$$gv | awk -F. '{if($$1>0 && $$2>=4 && $$3>=4) print "-DSSE4=1 -msse4.2"; else print "-DSSE4=0"}'))
endif
endif
