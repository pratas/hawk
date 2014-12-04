BIN    = .
CC     = gcc
CPLP   = -g -pg
#-----------------------------------------------------------------------------
CFLAGS = -O2 -Wall $(CPLP) -DPROGRESS #-DDEBUG
#-----------------------------------------------------------------------------
LIBS   = -lm
DEPS   = defs.h
PROGS  = $(BIN)/Hawk
OBJS   = mem.o misc.o args.o param.o dna.o hash.o bloom.o classes.o models.o \
         info.o bitio.o arith.o arith_aux.o
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
dna.o: dna.c dna.h $(DEPS)
	$(CC) -c $(CFLAGS) dna.c
hash.o: hash.c hash.h $(DEPS)
	$(CC) -c $(CFLAGS) hash.c
bloom.o: bloom.c bloom.h $(DEPS)
	$(CC) -c $(CFLAGS) bloom.c
classes.o: classes.c classes.h $(DEPS)
	$(CC) -c $(CFLAGS) classes.c
models.o: models.c models.h $(DEPS)
	$(CC) -c $(CFLAGS) models.c
info.o: info.c info.h $(DEPS)
	$(CC) -c $(CFLAGS) info.c
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
