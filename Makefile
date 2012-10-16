GCC44 := $(shell expr `gcc -dumpversion | cut -f1 -d.` \>= 4)$(shell expr `gcc -dumpversion | cut -f2 -d.` \>= 4)
DEBUG := 0
PROFILE := 0

ALL: mrsfast

LDFLAGS=-static
LIBS=-lz -lm -pthread -lpthread
CFLAGS=-O2 


ifeq "$(GCC44)" "11"
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

mrsfast: baseFAST.o MrsFAST.o Common.o CommandLineParser.o RefGenome.o HashTable.o Reads.o Output.o SnipReader.o HELP.o 
	gcc $^ -o $@ ${LDFLAGS} ${LIBS}
	rm -rf *.o

clean:
	rm -rf *.o

HELP.o:
	ld -r -b binary -o HELP.o HELP
