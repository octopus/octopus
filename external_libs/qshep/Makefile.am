# $Id$

noinst_LIBRARIES = libqshep.a

libqshep_a_SOURCES = qshep2d.f90 qshep3d.f90

# how to compile a f90 file :)
.f90.o:
	@FC@ @FCFLAGS@ -c @FCFLAGS_f90@ -o $@ $<


