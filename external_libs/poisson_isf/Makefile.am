# $Id: Makefile.am 1809 2006-02-09 20:39:38Z appel $

noinst_LIBRARIES = libpoisson_isf.a


libpoisson_isf_a_SOURCES = \
	Build_Kernel.f90  \
	fft3d.f90  \
	PSolver_Kernel.f90 \
	scaling_function.f90 \
	smooth.f90  \
        gequad.f90

# These two files are the parallel code, not used right now
#	pfft3d.f90 \
#	pconv.f90  

# how to compile a f90 file :)
.f90.o:
	@FC@ @FCFLAGS@ -c @FCFLAGS_f90@ -o $@ $<

scaling_function.o: lazy_100.inc  lazy_14.inc  \
	lazy_16.inc  lazy_20.inc  lazy_24.inc  \
	lazy_30.inc  lazy_40.inc  lazy_50.inc  \
	lazy_60.inc  lazy_8.inc
