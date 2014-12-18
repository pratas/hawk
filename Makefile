BIN    = .
CC     = gcc
CPLP   = #-g -pg
#-----------------------------------------------------------------------------
CFLAGS = -O3 -Wall $(CPLP) -DPROGRESS -DMEMORY #-DDEBUG -DREVERSE
#-----------------------------------------------------------------------------
LIBS   = -lm
DEPS   = defs.h
PROGS  = $(BIN)/Hawk
OBJS   = mem.o misc.o args.o param.o reads.o dna.o hash.o cch.o phash.o \
         classes.o sfcm.o models.o info.o gun.o bitio.o arith.o arith_aux.o
#-----------------------------------------------------------------------------
all:
	$(MAKE) progs
progs: $(PROGS)
$(BIN)/Hawk: hawk.c $(DEPS) $(OBJS)
	$(CC) $(CFLAGS) -o $(BIN)/Hawk hawk.c $(DEPS) $(OBJS) $(LIBS)
mem.o: mem.c mem.h $(DEPS)
	$(CC) -c $(CFLAGS) mem.c
misc.o: misc.c misc.h $(DEPS)
	$(CC) -c $(CFLAGS) misc.c
args.o: args.c args.h $(DEPS)
	$(CC) -c $(CFLAGS) args.c
param.o: param.c param.h $(DEPS)
	$(CC) -c $(CFLAGS) param.c
reads.o: reads.c reads.h $(DEPS)
	$(CC) -c $(CFLAGS) reads.c
dna.o: dna.c dna.h $(DEPS)
	$(CC) -c $(CFLAGS) dna.c
hash.o: hash.c hash.h $(DEPS)
	$(CC) -c $(CFLAGS) hash.c
phash.o: phash.c phash.h $(DEPS)
	$(CC) -c $(CFLAGS) phash.c
cch.o: cch.c cch.h $(DEPS)
	$(CC) -c $(CFLAGS) cch.c
classes.o: classes.c classes.h $(DEPS)
	$(CC) -c $(CFLAGS) classes.c
sfcm.o: sfcm.c sfcm.h $(DEPS)
	$(CC) -c $(CFLAGS) sfcm.c
models.o: models.c models.h $(DEPS)
	$(CC) -c $(CFLAGS) models.c
info.o: info.c info.h $(DEPS)
	$(CC) -c $(CFLAGS) info.c
gun.o: gun.c gun.h $(DEPS)
	$(CC) -c $(CFLAGS) gun.c
bitio.o: bitio.c bitio.h $(DEPS)
	$(CC) -c $(CFLAGS) bitio.c
arith.o: arith.c arith.h $(DEPS)
	$(CC) -c $(CFLAGS) arith.c
arith_aux.o: arith_aux.c arith_aux.h $(DEPS)
	$(CC) -c $(CFLAGS) arith_aux.c

clean:
	/bin/rm -f *.o
cleanall:
	/bin/rm -f *.o $(PROGS)
#-----------------------------------------------------------------------------
