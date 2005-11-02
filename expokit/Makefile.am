# $Id$

noinst_LIBRARIES = libexpokit.a

##libexpokit_a_SOURCES = expokit.f

##.f.o:
##	@FC@ @FCFLAGS@ -c -o $@ $<

libexpokit_a_SOURCES = expokit.f90

# how to compile a f90 file :)
.f90.o:
	@FC@ @FCFLAGS@ -c @FCFLAGS_f90@ -o $@ $<


