CC=mpicc
TEMP=.tmp
OBJS=nbc.o \
		 nbc_op.o \
		 dict.o \
		 hb_tree.o \
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
		 nbc_ialltoall.o \
		 nbc_fortran.o

OBJSF=$(shell for i in $(OBJS); do echo $$i | sed -e 's/\./_f\./'; done;)
#MELLANOX_LIBS=-L/usr/local/ibgd/driver/infinihost/lib/ -lvapi -lpthread -lmtl_common -lmosal -lmpga

#CFLAGS=-g -O3 -W -Wall -Winline -ansi -pedantic 
#CFLAGS=-g -Wall
CFLAGS=-O3 -Wall -I/usr/local/ibgd/driver/infinihost/include/
LDFLAGS=$(MELLANOX_LIBS) -lm

force: all

all: $(OBJS) libnbc.a 

%.o: %.c nbc.h 
	$(CC) $(CFLAGS)  -c $< -o $@ 

libnbc.a: $(OBJS)
	ar -rc $@ $(OBJS)

# prevent automatic deletion ...
.SECONDARY: 
	*m4

dist:
	for i in $$(ls -1 *c *h); do \
		cat COPYRIGHT $$i > $(TEMP); \
		mv $(TEMP) $$i; \
	done;	
	rm Makefile.orig
	mv README.dist README
	rm -rf CMakeFiles CMakeCache.txt cmake_install.cmake libdict-0.2.1.tar.gz
	rm -rf ib_main.c ib.c ib.h main*

clean: 
	rm -rf $(OBJS) nbc_test libnbc.a main.o nbc_f.h ib_test
	rm -rf $(OBJSF) *_f.c *_f.c.m4 libnbcf.a
