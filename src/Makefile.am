bin_PROGRAMS = @PROGNAME@
EXTRA_PROGRAMS = octopus octopus_1d

octopus_SOURCES = global.F90 liboct.F90 io.F90 pulpo.F90 units.F90 math.F90 \
	spline.F90 fft.F90 ps_xc.F90 ps_kb.F90 ps.F90 specie.F90 atom.F90 \
	mesh.F90 hartree.F90 states.F90 xc.F90 systm.F90 lasers.F90 hamiltonian.F90 \
	eigen.F90 mix.F90 lcao.F90 scf.F90 unocc.F90 static_pol.F90 \
	td.F90 run.F90 main.F90

noinst_HEADERS = config_F90.h \
	real.F90 complex.F90 undef.F90 \
	mesh1D.F90 mesh1D_create.F90 mesh1D_inc.F90 ps1D.F90 hartree1D.F90 \
	mesh3D.F90 mesh3D_create.F90 mesh3D_inc.F90 ps3D.F90 hartree3D.F90 \
	xc_pot.F90 xc_LDA.F90 xc_GGA.F90 xc_KLI.F90 xc_KLI_x.F90 xc_KLI_SIC.F90 \
	hamiltonian_inc.F90 h_external_pot.F90 h_forces.F90 \
	eigen_cg1.F90 eigen_cg2.F90 \
	states_inc.F90 \
	td_init.F90 td_rti.F90 \
	mix_inc.F90

octopus_LDADD = $(LIB_FFTW) $(LIB_GSL) $(LIB_LAPACK) ../liboct/liboct.a

octopus_1d_SOURCES = $(octopus_SOURCES)
octopus_1d_LDADD   = $(octopus_LDADD)

CLEANFILES = *~ *.bak *.mod config_F90.h
DISTCLEANFILES = Makefile.in

.F90.o:
	/lib/cpp -C $< > $*.f90
	$(F90) $(F90FLAGS)	-c -o $@ $*.f90
	@rm -f $*.f90

global.F90: config_F90.h
config_F90.h: config.h
	sed 's%/\*%!%g; s%\*/%%g' config.h >config_F90.h
../liboct/liboct.a:
	(cd ../liboct; $(MAKE))