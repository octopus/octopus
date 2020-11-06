# To be installed:
dynlibdir          = $(libdir)
noinst_LTLIBRARIES = $(flib_noinst)
lib_LTLIBRARIES    = $(flib_lib)
#module_DATA      = $(modules)
dynlib_DATA        = $(flib_dynlib)
noinst_DATA        = $(flib_linklib)

flib_noinst = libflib-1.la

# Source code.
#if HAVE_FC_FLUSH
#flush_files = flush.f90
#flush_object= flush.o
#else
flush_files = flush_fake.f90
flush_object= flush_fake.lo
#endif
if USE_MPI
mpi_source =
mpi_include =
mpi_inc=
else
mpi_source = MPIfake.f90
mpi_include = mpif.h
mpi_inc = -I$(builddir) 
endif
#if HAVE_FC_GET_COMMAND_ARGUMENT
getarg_files =
#else
#getarg_files = get_command_argument_fake.f90
#endif
libflib_1_la_SOURCES = utils.c utilsadd.c \
	yaml_strings.f90 dictionaries_base.f90 dictionaries.f90 \
	callbacks.f90 f_malloc.f90 \
	f_utils.f90 \
	yaml_output.f90 mem_profiling.f90 dynamic_memory.f90 \
	razero.f90 getadd.f90 time.f90 f_lib_highlev.f90 \
	yaml_parser.c yaml_parse.f90 \
	$(flush_files) $(mpi_source) $(getarg_files)

DYN_INCLUDES =  deallocate-inc.f90 allocate-inc.f90 \
	deallocate-c-inc.f90 allocate-c-inc.f90 \
	deallocate-profile-inc.f90 allocate-profile-inc.f90 \
	malloc_templates-inc.f90 \
	deallocate-multiple-inc.f90 deallocate-multiple-inc-ptr.f90 \
	f_memcpy-inc.f90 f_memcpy-base-inc.f90

F_MALLOC_INCLUDES = f_malloc-inc.f90 f_malloc-extra-inc.f90 \
	f_malloc-base-inc.f90 f_malloc-null-inc.f90 \
	f_malloc-total-inc.f90 f_malloc-bound-inc.f90 \
	f_malloc-bounds-inc.f90 f_malloc-simple-inc.f90 \
	f_malloc-ptr-inc.f90 f_malloc-check-inc.f90

YAML_INCLUDES = yaml_map-inc.f90 yaml_map-arr-inc.f90 \
	yaml_map-mat-inc.f90 yaml_open-inc.f90
STRINGS_INCLUDES = yaml_toa-inc.f90 yaml_toa-arr-inc.f90
ERR_INCLUDES = get_err-inc.f90
DICT_INCLUDES = error_handling.f90 $(ERR_INCLUDES) dict_add-inc.f90 dict_getvec-inc.f90 \
	set_arr-inc.f90 dict_cont-inc.f90 dict_cont_arr-inc.f90
GETADD_INCLUDES = getadd-c-inc.f90
NOOMP_INCLUDES = remove_omp-inc.f90 halt_omp-inc.f90

EXTRA_DIST = $(DYN_INCLUDES) $(YAML_INCLUDES) $(STRINGS_INCLUDES) $(ERR_INCLUDES) \
	$(DICT_INCLUDES) $(GETADD_INCLUDES) $(F_MALLOC_INCLUDES) $(NOOMP_INCLUDES) f_utils.inc

# Module files.
# if CAPITALIZE
# modules = YAML_STRINGS.@MODULE_EXT@ \
# 	DICTIONARIES.@MODULE_EXT@ \
# 	DICTIONARIES_BASE.@MODULE_EXT@ \
# 	YAML_OUTPUT.@MODULE_EXT@ \
# 	YAML_PARSE.@MODULE_EXT@ \
# 	METADATA_INTERFACES.@MODULE_EXT@ \
# 	MODULE_F_MALLOC.@MODULE_EXT@ \
# 	DYNAMIC_MEMORY.@MODULE_EXT@ \
# 	RANDOMDATA.@MODULE_EXT@ \
# 	MEMORY_PROFILING.@MODULE_EXT@ \
# 	F_UTILS.@MODULE_EXT@ \
# 	EXCEPTION_CALLBACKS.@MODULE_EXT@ \
# 	TIME_PROFILING.@MODULE_EXT@
# else
modules = yaml_strings.mod \
	dictionaries.mod \
	dictionaries_base.mod \
	yaml_output.mod \
	yaml_parse.mod \
	metadata_interfaces.mod \
	module_f_malloc.mod \
	dynamic_memory.mod \
	randomdata.mod \
	memory_profiling.mod \
	f_utils.mod \
	exception_callbacks.mod \
	time_profiling.mod
#endif

CLEANFILES = $(modules) *.mod $(mpi_include) memstatus.yaml

#AM_FCFLAGS = -I$(srcdir) $(mpi_inc) @LIBABINIT_INCLUDE@ @MPI_INCLUDE@
AM_FCFLAGS = -I. 
#-I$(srcdir) $(mpi_inc) 
#@MPI_INCLUDE@
AM_CFLAGS = @LIB_YAML_CFLAGS@

# Dependencies
$(modules): $(libflib_1_la_OBJECTS)
yaml_strings.lo: $(STRINGS_INCLUDES)
dictionaries_base.lo: yaml_strings.lo
dictionaries.lo: yaml_strings.lo dictionaries_base.lo $(DICT_INCLUDES) callbacks.lo
dictionaries_base.lo: yaml_strings.lo
callbacks.lo: $(NOOMP_INCLUDES)
#mem_profiling.lo : error_handling.lo
yaml_output.lo : dictionaries.lo f_utils.lo $(YAML_INCLUDES)
yaml_parse.lo : dictionaries.lo yaml_output.lo
#error_handling.lo: dictionaries.lo callbacks.lo yaml_output.lo $(ERR_INCLUDES)
f_malloc.lo: dictionaries.lo $(F_MALLOC_INCLUDES)
dynamic_memory.lo: getadd.lo yaml_strings.lo yaml_output.lo yaml_parse.lo mem_profiling.lo f_malloc.lo $(DYN_INCLUDES)
mem_profiling.lo: yaml_output.lo $(mpi_include) yaml_parse.lo
f_lib_highlev.lo: dictionaries.lo yaml_output.lo dynamic_memory.lo yaml_parse.lo time.lo
f_utils.lo: dictionaries.lo yaml_strings.lo f_utils.inc
getadd.lo: yaml_strings.lo $(GETADD_INCLUDES)
$(flush_object): yaml_output.lo dictionaries.lo

time.lo: $(NOOMP_INCLUDES) dynamic_memory.lo $(mpi_include) #this is temporary, mpif.h should only be in wrappers

#Fake mpif.h for serial compilation
mpif.h:
	touch mpif.h &&\
	echo "integer, parameter :: MPI_COMM_NULL=2,MPI_SUCCESS=0,MPI_SUM=1, MPI_COMM_WORLD=1" >> mpif.h &&\
	echo "integer, parameter :: MPI_DOUBLE_PRECISION=1, MPI_REAL=1, MPI_INTEGER=1" >> mpif.h &&\
	echo "integer, parameter :: MPI_STATUSES_IGNORE=1, MPI_LOGICAL=1" >> mpif.h &&\
	echo "integer, parameter :: MPI_MIN=1, MPI_MAX=1, MPI_CHARACTER=1, MPI_REAL8=1" >> mpif.h &&\
	echo "integer, parameter :: MPI_MAX_PROCESSOR_NAME=10, MPI_STATUS_SIZE=1,MPI_LAND=1" >> mpif.h &&\
	echo "integer, parameter :: MPI_REQUEST_NULL=1,MPI_STATUS_IGNORE=1" >> mpif.h &&\
	echo "integer, parameter :: mpi_tag_ub=1,mpi_address_kind=8,mpi_info_null=0" >> mpif.h &&\
	echo "integer, parameter :: mpi_mode_noprecede=0" >> mpif.h &&\
	echo "real(kind=8), external :: mpi_wtime" >> mpif.h


# Copy modules in a common include directory.
all: all-am 

