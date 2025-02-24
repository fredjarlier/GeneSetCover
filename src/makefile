CC=			/bioinfo/local/build/openmpi-1.8.3/bin/mpicc
CXX=		mpic++
CFLAGS=		-DBUG -g -Wall -o -O0
WRAP_MALLOC=-DUSE_MALLOC_WRAPPERS
AR=			ar
DFLAGS=		-DMPI -DHAVE_PTHREAD -D_FILE_OFFSET_BITS=64 -D_USE_FILE_OFFSET64 -D_LARGEFILE64_SOURCE -D_USE_MISC $(WRAP_MALLOC)
LOBJS=		
AOBJS=		 mergeSort.o mpi_globals.o parser.o tokenizer.o write.o genes.o duplicata.o biSort.o tree.o clusterCommunication.o read_communication.o condition2.o exon.o
			
PROG=		psort
INCLUDES=	
LIBS=		-lm -lz -lpthread -ltorque -lmpi -lopen-rte -lopen-pal -ldl -lnsl -lutil -ldl
SUBDIRS=	.

.SUFFIXES:.c .o .cc

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

psort:$(AOBJS) parallelMergeSort.o
		$(CC) $(CFLAGS) $(DFLAGS) $(AOBJS) parallelMergeSort.o -o $@ -L. $(LIBS)

#pbwa:libbwa.a $(AOBJS) main_serial_version.o
#		$(CC) $(CFLAGS) $(DFLAGS) $(AOBJS) main_serial_version.o -o $@ -L. -lbwa $(LIBS)

#bwamem-lite:libbwa.a example.o
#		$(CC) $(CFLAGS) $(DFLAGS) example.o -o $@ -L. -lbwa $(LIBS)

clean:
		rm -f gmon.out *.o a.out $(PROG) *~ *.a

depend:
	( LC_ALL=C ; export LC_ALL; makedepend -Y -- $(CFLAGS) $(DFLAGS) -- *.c )

# DO NOT DELETE THIS LINE -- make depend depends on it.

exon.o: exon.h
condition2.o: condition2.h
duplicata.o: duplicata.h
genes.o: genes.h
parser.o: parser.h tokenizer.h
biSort.o: parser.h mergeSort.h
main_parallel_version.o: parser.h 
write.o: write.h