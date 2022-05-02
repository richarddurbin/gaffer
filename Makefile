# makefile for gaffer developed on Richard's Mac

#CFLAGS= -O3
CFLAGS= -g -target arm64-apple-macos11				# for debugging
#CFLAGS= -03 -DOMP -fopenmp		# for OMP parallelisation - doesn't compile on Mac

ALL=seqconvert composition gaffer

DESTDIR=~/bin

all: $(ALL)

install:
	cp $(ALL) $(DESTDIR)

clean:
	$(RM) *.o *~ $(ALL)
	\rm -r *.dSYM

### object files

UTILS_OBJS=hash.o dict.o array.o utils.o
UTILS_HEADERS=utils.h array.h dict.h hash.h
$(UTILS_OBJS): utils.h $(UTILS_HEADERS)

SEQIO_OPTS = -DONEIO
SEQIO_LIBS = -lm -lz 

seqio.o: seqio.c seqio.h 
	$(CC) $(CFLAGS) $(SEQIO_OPTS) -c $^

ONElib.o: ONElib.c ONElib.h 
	$(CC) $(CFLAGS) -c $^

### programs

seqconvert: seqconvert.c seqio.o ONElib.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(SEQIO_LIBS)

composition: composition.c seqio.o ONElib.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(SEQIO_LIBS)

gaffer: gaffer.c seqio.o ONElib.o $(UTILS_OBJS)
	$(CC) $(CFLAGS) $^ -o $@ $(SEQIO_LIBS)

### end of file
