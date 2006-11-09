# $Id: Makefile.am 1809 2006-02-09 20:39:38Z appel $

noinst_LIBRARIES = libpoisson_isf.a


libpoisson_isf_a_SOURCES = \
	Build_Kernel.f90  \
	fft3d.f90  \
	pconv.f90  \
	pfft3d.f90  \
	PSolver_Kernel.f90 \
	scaling_function.f90 \
	smooth.f90  \
	dcopy.f  \
        gequad.f

# how to compile a f90 file :)
.f90.o:
	@FC@ @FCFLAGS@ -c @FCFLAGS_f90@ -o $@ $<


