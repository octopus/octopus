# To be installed:
noinst_LIBRARIES = libwrappers.a
#module_DATA = $(modules)

MPI_INCLUDES = allreduce-inc.f90 allreduce-arr-inc.f90 gather-inc.f90 \
	gather-inner-inc.f90 maxdiff-arr-inc.f90 maxdiff-inc.f90 \
	bcast-inc.f90 bcast-decl-inc.f90 bcast-decl-arr-inc.f90 \
	maxdiff-decl-inc.f90 maxdiff-end-inc.f90

if USE_MPI
#mpi_source =
mpi_include = $(MPI_INCLUDES)
else
#mpi_source = MPIfake.f90
mpi_include = $(MPI_INCLUDES) mpif.h
endif

EXTRA_DIST = $(MPI_INCLUDES)

#Fake mpif.h for serial compilation
mpif.h:
	touch mpif.h &&\
	echo "integer, parameter :: MPI_COMM_NULL=2, MPI_SUCCESS=0, MPI_SUM=1, MPI_COMM_WORLD=1" >> mpif.h &&\
	echo "integer, parameter :: MPI_DOUBLE_PRECISION=1, MPI_REAL=1, MPI_INTEGER=1, MPI_INTEGER8=1" >> mpif.h &&\
	echo "integer, parameter :: MPI_STATUSES_IGNORE=1, MPI_LOGICAL=1" >> mpif.h &&\
	echo "integer, parameter :: MPI_MIN=1, MPI_MAX=1, MPI_CHARACTER=1, MPI_REAL8=1" >> mpif.h &&\
	echo "integer, parameter :: MPI_MAX_PROCESSOR_NAME=10, MPI_STATUS_SIZE=1, MPI_LAND=1" >> mpif.h &&\
	echo "integer, parameter :: MPI_REQUEST_NULL=1, MPI_STATUS_IGNORE=1" >> mpif.h &&\
	echo "integer, parameter :: MPI_TAG_UB=1, MPI_ADDRESS_KIND=8, MPI_INFO_NULL=0" >> mpif.h &&\
	echo "integer, parameter :: MPI_MODE_NOPRECEDE=0" >> mpif.h &&\
	echo "real(kind=8), external :: mpi_wtime" >> mpif.h


# Source code.
#if USE_DGEMMSY
# DGEMMSY_SOURCES = interface_dgemmsy.f90 \
# 	dgemmsy/dgemmsy.c       dgemmsy/dgemmsy.h \
# 	dgemmsy/gemm_block_c.c  dgemmsy/gemm_block_c.h \
# 	dgemmsy/patterns.c      dgemmsy/patterns.h \
# 	dgemmsy/dgemmsy_utils.c dgemmsy/utils.h \
# 	dgemmsy/visitors.c      dgemmsy/visitors.h
#else
DGEMMSY_SOURCES = interface_dgemmsy_fake.f90
#endif
# if USE_CUDA_GPU
# CUBLAS_SOURCES = fortran.c
# else
CUBLAS_SOURCES = cublas_fake.f90 interface_cuda_fake.f90
#endif
#if USE_OCL
#OCL_SOURCES =
#else
OCL_SOURCES = interface_ocl_fake.f90
#endif
libwrappers_a_SOURCES = $(DGEMMSY_SOURCES) $(CUBLAS_SOURCES) $(OCL_SOURCES) \
	fft2d.f90 \
	fft3d.f90 \
	linalg.f90 \
	mpi.F90

# Dependencies
$(modules): $(libwrappers_a_OBJECTS)
interface_cuda_fake.o: linalg.o
interface_dgemmsy_fake.o: linalg.o

interface_dgemmsy_fake.o: linalg.o
mpi.o: $(mpi_include)


AM_FCFLAGS = -I$(top_builddir)/includes 
#@MPI_INCLUDE@
#AM_CPPFLAGS = @DGEMMSY_CPPFLAGS@
#if USE_CUDA_GPU
#AM_CPPFLAGS += -I@CUDA_PATH@/include @LIBSGPU_INCLUDE@
#endif

# Module files.
# if CAPITALIZE
# modules = MODULE_FFT_SG.mod \
# 	WRAPPER_LINALG.mod \
# 	WRAPPER_MPI.mod
# else
modules = module_fft_sg.mod \
	wrapper_linalg.mod \
	wrapper_mpi.mod
#endif
CLEANFILES = $(modules) mpif.h

# Copy modules in a common include directory.
all: all-am insmod
check: check-am insmod
insmod: $(modules)
	test -e "$(top_builddir)/includes" || $(INSTALL) -d -m 755 $(top_builddir)/includes
	for mod in "$(modules)" ; do \
	  $(INSTALL) -m 644  $$mod $(top_builddir)/includes ; \
	done
