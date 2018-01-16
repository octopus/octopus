## Process this file with automake to produce Makefile.in

OTHERS = bench.csh \
	acc_F.20-100.ref \
	acc_S-128.20-100.ref \
	acc_P.20-100.ref \
	README \
	INSTALL \
	$(REFERENCES) \
	$(EXERCISE)

REFERENCES= PS_Check.ref.yaml GPS_PCG.ref.yaml  PS_Basics.ref.yaml  PS_StressCheck.ref.yaml GPS_mPB.ref.yaml  GPS_PI.ref.yaml

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

EXTRA_DIST =  tols.yaml $(OTHERS)

check_PROGRAMS = PS_Check PS_Program PS_Exercise Generalized_PSolver Fock PS_StressCheck

AM_FCFLAGS = -I. -I$(top_srcdir)/src -I$(top_builddir)/src @LIB_FUTILE_CFLAGS@ @LIBETSFIO_INCLUDE@ @MPI_INCLUDE@
AM_LDFLAGS = -L$(top_builddir)/src
if BUILD_DYNAMIC_LIBS
AM_LDFLAGS += -Wl,-rpath=$(abs_top_builddir)/src -Wl,-rpath=$(DESTDIR)$(libdir)
psolver_library = $(top_builddir)/src/libPSolver-1.so.@PSOLVER_MINOR_VERSION@
else
psolver_library = $(top_builddir)/src/libPSolver-1.a @LIB_PSOLVER_DEPS@
endif

PS_Check_SOURCES = PS_Check.f90 functions.f90
PS_Check_LDADD = $(psolver_library) @LINALG_LIBS@

PS_Program_SOURCES = PS_Program.f90 functions.f90
PS_Program_LDADD = $(psolver_library) @LINALG_LIBS@

PS_StressCheck_SOURCES = PS_StressCheck.f90
PS_StressCheck_LDADD = $(psolver_library) @LINALG_LIBS@

PS_Exercise_SOURCES = PS_Exercise.f90
PS_Exercise_LDADD = $(psolver_library) @LINALG_LIBS@

Generalized_PSolver_SOURCES = Generalized_PSolver.f90
Generalized_PSolver_LDADD = $(psolver_library) @LINALG_LIBS@

Fock_SOURCES = Fock.f90 functions.f90
Fock_LDADD = $(psolver_library) @LINALG_LIBS@

if BUILD_DYNAMIC_LIBS
PS_Check_LDADD += @LIB_FUTILE_LIBS@
PS_Program_LDADD += @LIB_FUTILE_LIBS@
PS_StressCheck_LDADD += @LIB_FUTILE_LIBS@
PS_Exercise_LDADD += @LIB_FUTILE_LIBS@
Generalized_PSolver_LDADD += @LIB_FUTILE_LIBS@
Fock_LDADD += @LIB_FUTILE_LIBS@
endif

CLEANFILES = malloc.prc time.prc *.out.yaml *.report.yaml report_remarks *dat fort* log.yaml time.yaml

#PS_Check.out.yaml: PS_Check Generalized_PSolver
#	-$(run_parallel) ./PS_Check -n [57,48,63] -g F > PS_Check.out.yaml
#	-$(run_parallel) ./PS_Check -n [64,64,64] -g P >> PS_Check.out.yaml
#	-$(run_parallel) ./PS_Check -n [32,64,48] -g S >> PS_Check.out.yaml
#	-$(run_parallel) ./Generalized_PSolver -n [57,49,63] -m PI -i '{environment: {input_guess: No}}'
#	cat log.yaml >> PS_Check.out.yaml
#	-$(run_parallel) ./Generalized_PSolver -n [57,49,63] -m PCG -i '{environment: {input_guess: No}}'
#	cat log.yaml >> PS_Check.out.yaml
#	cp PS_Check.out.yaml tmp
#	cat tmp | grep -v 'Unable to read mpd.hosts' > PS_Check.out.yaml
#	rm tmp

check: regtest report #PS_Basics.out.yaml PS_Check.out.yaml PS_Check.ref.yaml tols.yaml report

report: PS_Basics.report.yaml PS_Check.report.yaml PS_Basics.report.yaml
	@if test $(MAKELEVEL) = 0 ; then python $(pythondir)/report.py ; fi

PS_Basics.report.yaml: PS_Basics.ref.yaml PS_Basics.out.yaml
	python $(pythondir)/fldiff_yaml.py \
		-r $< -d PS_Basics.out.yaml -t $(srcdir)/tols.yaml \
	        --output $@

PS_Check.report.yaml: PS_Check.ref.yaml PS_Check.out.yaml
	python $(pythondir)/fldiff_yaml.py \
		-r $< -d PS_Check.out.yaml -t $(srcdir)/tols.yaml \
	        --output $@

regtest: PS_Program PS_Exercise PS_Check Generalized_PSolver PS_StressCheck
	python $(pythondir)/f_regtest.py --fldiff=$(pythondir)/fldiff_yaml.py \
	       --tols=$(srcdir)/tols.yaml --srcdir=$(srcdir) \
	       --exclusive=${CHECK_ONLY}


define F_REGTEST_INSTRUCTIONS
- PS_Basics:
   binary: PS_Program
   runs:
    # broken - $(run_parallel) ./PS_Program -n [57,48,63] -g F > PS_Basics.out.yaml
    - $(run_parallel) ./PS_Program -n [64,64,64] -g P
    - $(run_parallel) ./PS_Program -n [32,64,48] -g S
    - $(run_parallel) ./PS_Program -n [32,64,48] -g S -d [90,89,90]
    - $(run_parallel) ./PS_Program -n [32,64,48] -g S -d [90,60,90]
    # broken - $(run_parallel) ./PS_Program -n [32,64,48] -g W >> PS_Basics.out.yaml
- PS_Check:
   runs:
    - $(run_parallel) ./PS_Check -n [57,48,63] -g F
    - $(run_parallel) ./PS_Check -n [64,64,64] -g P
    - $(run_parallel) ./PS_Check -n [32,64,48] -g S
- GPS_PCG:
   binary: Generalized_PSolver
   runs:
    # PCG from psolver, ortho and check also some dsurfdrho subroutines with -c option.
    - $(run_parallel) ./Generalized_PSolver -n [57,49,63] -m PCG -c yes -l no
    # PCG from psolver, ortho with internal cavity vector built from eps in pkernel_set_epsilon.
    - $(run_parallel) ./Generalized_PSolver -n [57,49,63] -m PCG -b yes -l no
    # PCG from psolver, non-ortho
    - $(run_parallel) ./Generalized_PSolver -n [60,78,70] -m PCG -f [90,41,90] -g S -l no
    # PCG from Generalized_PSolver, non-ortho
    - $(run_parallel) ./Generalized_PSolver -n [60,78,70] -m VAC -e 5 -g S -f [90,41,90] -l no
- GPS_PI:
   binary: Generalized_PSolver
   runs:
    # PI from psolver, ortho
    - $(run_parallel) ./Generalized_PSolver -n [57,49,63] -m PI -l no
    # PI from psolver, ortho with internal cavity vector built from eps in pkernel_set_epsilon.
    - $(run_parallel) ./Generalized_PSolver -n [57,49,63] -m PI -b yes -l no
    # PI from Generalized_PSolver, non-ortho
    - $(run_parallel) ./Generalized_PSolver -n [60,78,70] -m VAC -e 6 -g S -f [90,41,90] -l no
- GPS_mPB:
   binary: Generalized_PSolver
   runs:
    # mPB from psolver, ortho
    - "$(run_parallel) ./Generalized_PSolver -n [60,78,70] -m PCG -g S -p yes -i '{environment: {minres: 1.e-6, pb_minres: 1.e-6}}' -l no"
    # mPB from Generalized_PSolver, non-ortho
    - $(run_parallel) ./Generalized_PSolver -n [60,78,70] -m VAC -e 8 -g S -f [90,41,90] -p yes -l no
- PS_StressCheck:
   runs:
    # Check stress on x,y,z direction individually
    - $(run_parallel) ./PS_StressCheck -n 50 -g P -s 50
    # Check stress on x,y,z direction concurrently
    - $(run_parallel) ./PS_StressCheck -n 50 -g P -s 50 -v yes
endef
export F_REGTEST_INSTRUCTIONS
