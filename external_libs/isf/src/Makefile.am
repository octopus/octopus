# To be installed:
dynlibdir        = $(libdir)
noinst_LIBRARIES = $(psolver_noinst)
lib_LIBRARIES    = $(psolver_lib)
# module_DATA      = $(modules)
dynlib_DATA      = $(psolver_dynlib)
noinst_DATA      = $(psolver_linklib)


# Static or dynamic library
psolver_dynamic_library = libisf.so.@BIGDFT_MINOR_VERSION@.0.@BIGDFT_MICRO_VERSION@
psolver_library_soname  = libisf.so.@BIGDFT_MINOR_VERSION@

# We install a PSolver library only if BigDFT is deactivated.
#if BUILD_LIB_BIGDFT
#psolver_lib    =
#psolver_dynlib =
  psolver_noinst = libisf.a
#else
#psolver_lib    = libPSolver-1.a
# if BUILD_DYNAMIC_LIBS
# psolver_dynlib = $(psolver_dynamic_library)
# psolver_linklib = $(psolver_library_soname)
# else
# psolver_dynlib =
# psolver_linklib =
# endif
# psolver_noinst =
# endif


# Code source.
# if USE_CUDA_GPU
# CUDA_SOURCES = cufft.cu
# else
CUDA_SOURCES =
#endif
libisf_a_SOURCES = Build_Kernel.f90 \
	scaling_function.f90 \
	PSolver_Base_new.f90 \
	wofz.f90 \
	$(CUDA_SOURCES) \
	Poisson_Solver.f90

# Modules.
#if CAPITALIZE
#modules = POISSON_SOLVER.@MODULE_EXT@
#else
modules = poisson_solver.mod
#endif

CLEANFILES = $(modules) $(psolver_dynamic_library) $(psolver_library_soname) libisf.so
SOURCES_MOD = PSolver_Main.f90 createKernel.f90
EXTRA_DIST = \
	gaussfit_wires.inc\
	gaussfit_Yukawa.inc\
	lazy_8.inc\
	lazy_14.inc\
	lazy_16.inc\
	lazy_20.inc\
	lazy_24.inc\
	lazy_30.inc\
	lazy_40.inc\
	lazy_50.inc\
	lazy_60.inc\
	lazy_100.inc\
	lazy_ISF_8_2048.inc\
	lazy_ISF_14_2048.inc\
	lazy_ISF_16_2048.preinc\
	lazy_ISF_16_2048.inc\
	lazy_ISF_20_2048.inc\
	lazy_ISF_24_2048.inc\
	lazy_ISF_30_2048.inc\
	lazy_ISF_40_2048.inc \
	lazy_ISF_50_2048.inc\
	lazy_ISF_60_2048.inc\
	lazy_ISF_100_2048.inc\
	perfdata.inc \
	configure.inc.in \
	$(SOURCES_MOD)

# Build targets.
# if USE_CUDA_GPU
# CONV_CHECK = conv_check_fft
# else
CONV_CHECK = conv_check_fft
# endif
EXTRA_PROGRAMS = $(CONV_CHECK)
conv_check_fft_SOURCES = conv_check_fft.f90
conv_check_fft_LDADD = libisf.a $(top_builddir)/wrappers/libwrappers.a $(top_builddir)/flib/src/libflib-1.a @LIBS_EXT_DEPS@
#conv_check_fftw_SOURCES = conv_check_fftw.f90 fftw.cpp
#conv_check_fftw_LDADD = @LIBS_EXT_DEPS@ -lfftw3 -lfftw3_threads

AM_FCFLAGS = -I$(top_builddir)/includes -I. -I$(srcdir) 
#@MPI_INCLUDE@
AM_CPPFLAGS = -I. -I$(top_builddir) 
#-I@CUDA_PATH@/include/ 
AM_LDFLAGS = -L$(top_builddir)/src
# if BUILD_DYNAMIC_LIBS
# AM_LDFLAGS += -Wl,-rpath=$(abs_top_builddir)/src
# else
LIBS += @LIBS_EXT_DEPS@
# endif

SUFFIXES = .cu

NVCC = @NVCC@
NVCC_FLAGS = -arch sm_13 @NVCC_FLAGS@

# Dynamic library building
$(psolver_dynamic_library): $(libisf_a_OBJECTS)
	$(FC) $(FCFLAGS) -shared $(LDFLAGS) -Wl,-soname=$(psolver_library_soname) -o $@ $^ $(top_builddir)/wrappers/libwrappers.a $(top_builddir)/flib/src/libflib-1.a @LIBS_EXT_DEPS@
	@chmod a+x $@

$(psolver_library_soname): $(psolver_dynamic_library)
	ln -fs $^ $@
	ln -fs $@ libisf.so

install-data-hook:
	if test -n "$(psolver_dynlib)" ; then \
	  cd $(DESTDIR)$(dynlibdir) ; \
	  chmod a+x $(psolver_dynamic_library) ; \
	  ln -fs $(psolver_dynamic_library) $(psolver_library_soname) ; \
	  ln -fs $(psolver_library_soname) libisf.so ; \
	fi

#dependencies
$(modules): $(libisf_a_OBJECTS)
Poisson_Solver.o: $(SOURCES_MOD) configure.inc
conv_check_fft.o: cufft.o
conv_check_fftw.o: fftw.o
wofz.o Build_Kernel.o scaling_function.o PSolver_Base_new.o: Poisson_Solver.o

.cu.o:
	$(NVCC) $(NVCC_FLAGS) $(AM_CPPFLAGS) -c $< -o $@

.cpp.o:
	icc -openmp -lfftw3 -lfftw3_threads -c $< -o $@

# Copy modules in a common include directory.
all: all-am insmod
check: check-am insmod
insmod: $(modules)
	test -e "$(top_builddir)/includes" || $(INSTALL) -d -m 755 $(top_builddir)/includes
	for mod in "$(modules)" ; do \
	  $(INSTALL) -m 644  $$mod $(top_builddir)/includes ; \
	done
