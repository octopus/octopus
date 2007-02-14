CC=mpicc
TEMP=.tmp
OBJS=nbc.o \
		 nbc_op.o \
		 nbc_ialltoallv.o \
		 nbc_ibarrier.o \
		 nbc_ibcast.o \
		 nbc_igather.o \
		 nbc_igatherv.o \
		 nbc_iscatter.o \
		 nbc_iscatterv.o \
		 nbc_iallgather.o \
		 nbc_iallgatherv.o \
		 nbc_ireduce.o \
		 nbc_ireduce_scatter.o \
		 nbc_iallreduce.o \
		 nbc_iscan.o \
		 nbc_ialltoall.o
OBJSF=$(shell for i in $(OBJS); do echo $$i | sed -e 's/\./_f\./'; done;)

#CFLAGS=-g -W -Wall -Winline -ansi -pedantic 
#CFLAGS=-g -Wall
CFLAGS=-O3

force: all

all: $(OBJS) nbc_test libnbc.a

%.o: %.c nbc.h
	$(CC) $(CFLAGS)  -c $< -o $@ 

nbc_test: libnbc.a main.c
	$(CC) $(CFLAGS) main.c libnbc.a -o nbc_test -lm

libnbc.a: $(OBJS)
	ar -rc $@ *.o

%_f.c.m4: %.c
	cat $< | sed -e 's/nbc.h/nbc_f.h/' > $(TEMP)
	cat c2f.m4 $(TEMP) > $@
	rm $(TEMP)

nbc_f.h: nbc.h
	cat c2f.m4 $< > $(TEMP)
	m4 $(TEMP) > $@

%_f.c: %_f.c.m4 nbc_f.h
	m4 $< > $@

libnbcf.a: $(OBJSF)
	ar -rc $@ *_f.o

# prevent automatic deletion ...
.SECONDARY: 
	*m4

dist:
	for i in $$(ls -1 *c); do \
		cat COPYRIGHT $$i > $(TEMP); \
		mv $(TEMP) $$i; \
	done;	
	mv $(TEMP) Makefile 
	mv README.dist README

clean: 
	rm -rf $(OBJS) nbc_test libnbc.a main.o nbc_f.h
	rm -rf $(OBJSF) *_f.c *_f.c.m4 libnbcf.a
