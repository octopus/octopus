# To be installed:
dynlibdir        = $(libdir)
lib_LIBRARIES    = libPSolver-1.a
module_DATA      = $(modules)
dynlib_DATA      = $(psolver_dynlib)
noinst_DATA      = $(psolver_linklib)

# Static or dynamic library
psolver_dynamic_library = libPSolver-1.so.@PSOLVER_MINOR_VERSION@.0.@PSOLVER_MICRO_VERSION@
psolver_library_soname  = libPSolver-1.so.@PSOLVER_MINOR_VERSION@

if BUILD_DYNAMIC_LIBS
psolver_dynlib = $(psolver_dynamic_library)
psolver_linklib = $(psolver_library_soname)
psolver_library = $(psolver_library_soname)
else
psolver_dynlib =
psolver_linklib =
psolver_library = libPSolver-1.a @LIB_PSOLVER_DEPS@
endif

# C objects
c_level = PS_input_variables_fetch.c

PS_input_variables_fetch.o: PS_input_variables_definition-inc.h
#build the fortran include file with the stream of input variables into a local string
#the sed command below will transform the entire yaml document into a c string, to be passed to fortran
PS_input_variables_definition-inc.h: PS_input_variables_definition.yaml
	$(SED) -e "s/^/\"/;s/$$/\\\n\"/" $< > $@

if HAVE_ETSF_IO
etsf_files = IOetsf.f90
etsf_dep   = IOetsf.o
else
etsf_files = IOetsf_fake.f90
etsf_dep   = IOetsf_fake.o
endif

fdder_file=FDder.f90

# Code source.
if USE_CUDA_GPU
CUDA_SOURCES = cufft.cu
else
CUDA_SOURCES = cufft_fake.f90
endif
libPSolver_1_a_SOURCES = Build_Kernel.f90 \
	scaling_function.f90 \
	PSolver_Base_new.f90 \
	wofz.f90 \
	$(CUDA_SOURCES) \
	Poisson_Solver.f90 \
	environment.f90 \
	PStypes.f90 \
	PSbox.f90 \
	PSbase.f90 \
	exctx_calculation.f90 \
	mp_quadrature.f90 \
	$(fdder_file) \
	IObox.f90 \
	$(c_level) \
	$(etsf_files)

# Modules.
if CAPITALIZE
modules = POISSON_SOLVER.@MODULE_EXT@ \
	PSOLVER_ENVIRONMENT.@MODULE_EXT@ \
	PSTYPES.@MODULE_EXT@ \
	PSBASE.@MODULE_EXT@ \
	PSBOX.@MODULE_EXT@ \
	FDDER.@MODULE_EXT@ \
	IOBOXETSF.@MODULE_EXT@ \
	MULTIPOLE_PRESERVING.@MODULE_EXT@ \
	IOBOX.@MODULE_EXT@
else
modules = poisson_solver.@MODULE_EXT@ \
	psolver_environment.@MODULE_EXT@ \
	pstypes.@MODULE_EXT@ \
	psbase.@MODULE_EXT@ \
	psbox.@MODULE_EXT@ \
	fdder.@MODULE_EXT@ \
	ioboxetsf.@MODULE_EXT@ \
	multipole_preserving.@MODULE_EXT@ \
	iobox.@MODULE_EXT@
endif

CLEANFILES = *.@MODULE_EXT@ $(psolver_dynamic_library) $(psolver_library_soname) libPSolver-1.so
SOURCES_MOD = PSolver_Main.f90 createKernel.f90
EXTRA_DIST = \
	gaussfit_wires.inc\
	gaussfit_Yukawa.inc\
	lazy_100.inc\
	lazy_14.inc\
	lazy_16.inc\
	lazy_20.inc\
	lazy_24.inc\
	lazy_2.inc\
	lazy_30.inc\
	lazy_40.inc\
	lazy_4.inc\
	lazy_50.inc\
	lazy_60.inc\
	lazy_6.inc\
	lazy_8.inc\
	lazy_ISF_100_2048.inc\
	lazy_ISF_14_2048.inc\
	lazy_ISF_16_2048.inc\
	lazy_ISF_16_2048.preinc\
	lazy_ISF_20_2048.inc\
	lazy_ISF_24_2048.inc\
	lazy_ISF_30_2048.inc\
	lazy_ISF_40_2048.inc\
	lazy_ISF_50_2048.inc\
	lazy_ISF_60_2048.inc\
	lazy_ISF_8_2048.inc\
	perfdata.inc\
	configure.inc.in\
	FiniteDiffCorff.inc\
	FiniteDiffCorff_2der.inc\
	PS_input_variables_definition.yaml \
	$(SOURCES_MOD)

# Build targets.
if USE_CUDA_GPU
CONV_CHECK = conv_check_fft
else
CONV_CHECK = conv_check_fftw
endif
EXTRA_PROGRAMS = $(CONV_CHECK)
conv_check_fft_SOURCES = conv_check_fft.f90
conv_check_fft_LDADD = $(psolver_library)
conv_check_fftw_SOURCES = conv_check_fftw.f90 fftw.cpp
conv_check_fftw_LDADD = $(psolver_library) -lfftw3 -lfftw3_threads

AM_FCFLAGS = -I. -I$(srcdir) @LIB_FUTILE_CFLAGS@ @LIBETSFIO_INCLUDE@ @MPI_INCLUDE@
AM_CPPFLAGS = -I. -I$(top_builddir) -I@CUDA_PATH@/include/ 
AM_LDFLAGS = -L$(top_builddir)/src

if BUILD_DYNAMIC_LIBS
LD_RPATH=-Wl,-rpath=$(DESTDIR)$(dynlibdir)
AM_LDFLAGS += $(LD_RPATH)
endif

SUFFIXES = .cu

NVCC = @NVCC@
NVCC_FLAGS = @NVCC_FLAGS@

# Dynamic library building
$(psolver_dynamic_library): $(libPSolver_1_a_OBJECTS)
	$(FC) $(FCFLAGS) -shared $(LDFLAGS) -Wl,-soname=$(psolver_library_soname) $(LD_RPATH) -o $@ $^ @LIB_PSOLVER_DEPS@ $(LIBS)
	@chmod a+x $@

$(psolver_library_soname): $(psolver_dynamic_library)
	ln -fs $^ $@
	ln -fs $@ libPSolver-1.so

install-data-hook:
	if test -n "$(psolver_dynlib)" ; then \
	  cd $(DESTDIR)$(dynlibdir) ; \
	  chmod a+x $(psolver_dynamic_library) ; \
	  ln -fs $(psolver_dynamic_library) $(psolver_library_soname) ; \
	  ln -fs $(psolver_library_soname) libPSolver-1.so ; \
	fi

#dependencies
$(modules): $(libPSolver_1_a_OBJECTS)
Poisson_Solver.o: $(SOURCES_MOD) configure.inc environment.o PStypes.o PSbox.o PSbase.o FDder.o
PStypes.o: environment.o FDder.o
PSbox.o: PStypes.o PSbase.o IObox.o
FDder.o: PSbase.o
IObox.o: PSbase.o $(etsf_dep)
environment.o: PSbase.o
conv_check_fft.o: cufft.o
conv_check_fftw.o: fftw.o
IOetsf.o: PSbase.o
IOetsf_fake.o: PSbase.o
exctx_calculation.o wofz.o Build_Kernel.o scaling_function.o PSolver_Base_new.o: Poisson_Solver.o

.cu.o:
	$(NVCC) $(NVCC_FLAGS) $(AM_CPPFLAGS) -c $< -o $@

.cpp.o:
	icc -openmp -lfftw3 -lfftw3_threads -c $< -o $@
