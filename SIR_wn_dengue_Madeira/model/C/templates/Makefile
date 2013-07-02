CC=gcc
CFLAGS=-fopenmp -std=gnu99 -O3 -DGSL_RANGE_CHECK_OFF
LDFLAGS=-lplom -lplomtpl -lm -lgsl -lgslcblas -ljansson -lzmq
LIB=libplomtpl.a
SRC= $(wildcard *.c)
OBJ= $(SRC:.c=.o)
LIBDIR=../lib

EXEC=simul simplex smc mif pmcmc kalman ksimplex kmcmc worker

all: $(LIB)

install: $(EXEC)
	mv $(EXEC) ../../

$(LIB): $(OBJ)
	ar -rcs $@ $(OBJ)

%.o: %.c
	$(CC) $(CFLAGS) -o $@ -c $<

simul: libplomtpl.a
	$(CC) $(CFLAGS) -L. -o $@ -lplom_simul $(LDFLAGS)

simplex: libplomtpl.a
	$(CC) $(CFLAGS) -L. -o $@ -lplom_simplex $(LDFLAGS)

smc: libplomtpl.a
	$(CC) $(CFLAGS) -L. -o $@ -lplom_smc $(LDFLAGS) 

mif: libplomtpl.a
	$(CC) $(CFLAGS) -L. -o $@ -lplom_mif $(LDFLAGS)  

pmcmc: libplomtpl.a
	$(CC) $(CFLAGS) -L. -o $@ -lplom_pmcmc $(LDFLAGS)  

kalman: libplomtpl.a
	$(CC) $(CFLAGS) -L. -o $@ -lplom_kalman $(LDFLAGS) 

ksimplex: libplomtpl.a
	$(CC) $(CFLAGS) -L. -o $@ -lplom_ksimplex $(LDFLAGS)

kmcmc: libplomtpl.a
	$(CC) $(CFLAGS) -L. -o $@ -lplom_kmcmc $(LDFLAGS)

worker: libplomtpl.a
	$(CC) $(CFLAGS) -L. -o $@ -lplom_worker $(LDFLAGS) 

.PHONY: clean

clean:
	rm *.o $(LIB)

uninstall: clean
	rm ../../{simul,simplex,smc,mif,pmcmc,kalman,ksimplex,kmcmc,worker}
