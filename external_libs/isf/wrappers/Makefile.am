# To be installed:
noinst_LIBRARIES = libwrappers.a
#module_DATA = $(modules)

MPI_INCLUDES = allreduce-inc.f90 allreduce-arr-inc.f90 gather-inc.f90 \
	gather-inner-inc.f90 maxdiff-arr-inc.f90 maxdiff-inc.f90 \
	bcast-inc.f90 bcast-decl-inc.f90 bcast-decl-arr-inc.f90 \
	maxdiff-decl-inc.f90 maxdiff-end-inc.f90

if USE_MPI
  mpi_include = $(MPI_INCLUDES)
else
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

DGEMMSY_SOURCES = interface_dgemmsy_fake.f90
CUBLAS_SOURCES = cublas_fake.f90 interface_cuda_fake.f90
OCL_SOURCES = interface_ocl_fake.f90

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

# Compilation is a two-step process: first we preprocess F90 files
# to generate _oct.f90 files. Then, we compile this _oct.f90 into
# an object file and delete the intermediate file.
.F90.o:
	@FCCPP@ @CPPFLAGS@ $(DEFS) $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) -I. $< > $*_oct.f90
	$(top_srcdir)/build/preprocess.pl $*_oct.f90 "@DEBUG@" "@F90_ACCEPTS_LINE_NUMBERS@" "@F90_FORALL@"
	@FC@ @FCFLAGS@ $(AM_FCFLAGS) -c @FCFLAGS_f90@ -o $@ $*_oct.f90
	@rm -f $*_oct.f90

#AM_FCFLAGS = -I$(top_builddir)/includes 
AM_FCFLAGS = -I../flib -I.
modules = module_fft_sg.mod \
	wrapper_linalg.mod \
	wrapper_mpi.mod

CLEANFILES = $(modules) mpif.h

# Copy modules in a common include directory.
all: all-am 

