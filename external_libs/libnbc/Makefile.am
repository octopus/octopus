# $Id: Makefile.am 1809 2006-02-09 20:39:38Z appel $

if COMPILE_LIBNBC
  noinst_LIBRARIES = libnbc.a
endif

libnbc_a_SOURCES = \
	nbc.h \
	nbc.c \
	nbc_iallgather.c \
	nbc_iallgatherv.c \
	nbc_iallreduce.c \
	nbc_ialltoall.c \
	nbc_ialltoallv.c \
	nbc_ibarrier.c \
	nbc_ibcast.c \
	nbc_igather.c \
	nbc_igatherv.c \
	nbc_ireduce.c \
	nbc_ireduce_scatter.c \
	nbc_iscan.c \
	nbc_iscatter.c \
	nbc_iscatterv.c \
	nbc_op.c \
	nbcf.c

CLEANFILES = *~ *.bak
DISTCLEANFILES = Makefile.in
