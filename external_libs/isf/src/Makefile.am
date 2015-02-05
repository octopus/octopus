# To be installed:
dynlibdir        = $(libdir)
noinst_LIBRARIES = $(psolver_noinst)
lib_LIBRARIES    = $(psolver_lib)
# module_DATA      = $(modules)
dynlib_DATA      = $(psolver_dynlib)
noinst_DATA      = $(psolver_linklib)


psolver_noinst = libisf.a

libisf_a_SOURCES = Build_Kernel.f90 \
	scaling_function.f90 \
	PSolver_Base_new.f90 \
	wofz.f90 \
	Poisson_Solver.f90

modules = poisson_solver.mod

CLEANFILES = $(modules) $(psolver_dynamic_library) $(psolver_library_soname)
SOURCES_MOD = PSolver_Main.f90 createKernel.f90
noinst_HEADERS = $(SOURCES_MOD)

#AM_FCFLAGS = -I$(top_builddir)/includes -I. -I$(srcdir) 
AM_FCFLAGS = -I../wrappers -I../flib -I.  
AM_CPPFLAGS = -I. 
#-I$(top_builddir) 
AM_LDFLAGS = -L$(top_builddir)/src
LIBS += @LIBS_EXT_DEPS@

#dependencies
$(modules): $(libisf_a_OBJECTS)
Poisson_Solver.o: $(SOURCES_MOD) configure.inc
wofz.o Build_Kernel.o scaling_function.o PSolver_Base_new.o: Poisson_Solver.o
