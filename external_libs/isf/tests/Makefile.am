## Process this file with automake to produce Makefile.in

OTHERS = bench.csh \
	acc_F.20-100.ref \
	acc_S-128.20-100.ref \
	acc_P.20-100.ref \
	README \
	INSTALL \
	$(EXERCISE)

EXERCISE = Exercise/plot.gnuplot \
	Exercise/PS_Exercise.tex \
	Exercise/accF.pdf \
	Exercise/accS.pdf \
	Exercise/fandg.pdf \
	Exercise/PvsF.pdf \
	Exercise/rhoSpotS.pdf \
	Exercise/accP.pdf \
	Exercise/emepot.pdf \
	Exercise/rhoFpotF.pdf

EXTRA_DIST = PS_Check.ref.yaml tols.yaml $(OTHERS)

check_PROGRAMS = PS_Check PS_Program PS_Exercise

if BUILD_LIBYAML
LD_LIBRARY_PATH := ${LD_LIBRARY_PATH}:$(abs_top_builddir)/yaml-0.1.4/src/.libs
PYTHONPATH := ${PYTHONPATH}:`ls -d $(abs_top_builddir)/PyYAML-3.10/build/lib.*`
endif

AM_FCFLAGS = -I$(top_builddir)/includes -I$(top_srcdir)/PSolver/src -I. @MPI_INCLUDE@
AM_LDFLAGS = -L$(top_builddir)/PSolver/src

f_libraries = $(top_builddir)/wrappers/libwrappers.a $(top_builddir)/flib/src/libflib-1.a

if BUILD_LIB_BIGDFT
psolver_libraries = $(top_builddir)/PSolver/src/libPSolver-1.a $(f_libraries)
LIBS += @LIBS_EXT_DEPS@
else
if BUILD_DYNAMIC_LIBS
psolver_libraries = -lPSolver-1
AM_LDFLAGS += -Wl,-rpath=$(abs_top_builddir)/PSolver/src
else
psolver_libraries = $(top_builddir)/PSolver/src/libPSolver-1.a $(f_libraries)
LIBS += @LIBS_EXT_DEPS@
endif
endif

PS_Check_SOURCES = PS_Check.f90
PS_Check_LDADD = $(psolver_libraries)

PS_Program_SOURCES = PS_Program.f90
PS_Program_LDADD = $(psolver_libraries)

PS_Exercise_SOURCES = PS_Exercise.f90
PS_Exercise_LDADD = $(psolver_libraries)

CLEANFILES = malloc.prc time.prc PS_Check.out.yaml PSolver.report.yaml

PS_Check.out.yaml: PS_Check PS_Program PS_Exercise 
	-$(run_parallel) ./PS_Check 57 48 63  F  > PS_Check.out.yaml
	-$(run_parallel) ./PS_Check 64 64 64  P >> PS_Check.out.yaml
	-$(run_parallel) ./PS_Check 32 64 48  S >> PS_Check.out.yaml
	cp PS_Check.out.yaml tmp
	cat tmp | grep -v 'Unable to read mpd.hosts' > PS_Check.out.yaml
	rm tmp

check: PS_Check.out.yaml
	export PYTHONPATH=${PYTHONPATH} ; \
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH} ; \
	python $(top_srcdir)/tests/fldiff_yaml.py \
		-r $(abs_srcdir)/PS_Check.ref.yaml -d PS_Check.out.yaml -t $(abs_srcdir)/tols.yaml \
	        --output PSolver.report.yaml

#	python $(top_srcdir)/tests/fldiff.py --mode=psolver --discrepancy=2.e-8 \
#		PS_Check.out $(srcdir)/PS_Check.out.ref | tee fldiff.report
