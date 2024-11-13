# Common definitions
CC = gcc

# Compiler flags, paths and libraries
CFLAGS = -std=c11 -pedantic -Wall -Wextra -O3 -ggdb -gdwarf-4
LFLAGS = $(CFLAGS)
LIBS   = -lm -lpthread

TGTS = partdiff-posix partdiff-sep

# Targets ...
all: partdiff-posix partdiff-seq

partdiff-seq: partdiff-seq.o askparams.o Makefile
	$(CC) $(LFLAGS) -o $@ partdiff-seq.o askparams.o $(LIBS)

partdiff-posix: partdiff-posix.o askparams.o Makefile
	$(CC) $(LFLAGS) -o $@ partdiff-posix.o askparams.o $(LIBS)

partdiff-seq.o: partdiff.c Makefile
	$(CC) -c $(CFLAGS) -D SEQUENTIAL  partdiff.c -o partdiff-seq.o

partdiff-posix.o: partdiff.c Makefile
	$(CC) -c $(CFLAGS) partdiff.c -o partdiff-posix.o

askparams.o: askparams.c Makefile

# Rule to create *.o from *.c
%.o: %.c
	$(CC) -c $(CFLAGS) $*.c

clean:
	$(RM) $(OBJS)
	$(RM) $(TGTS)
