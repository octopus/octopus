bin_PROGRAMS = @PROGNAME@ strength-function hs-from-mult hs-from-acc \
	excitations broad function-to-dx phonon-spectrum nbo2xyz
EXTRA_PROGRAMS = octopus octopus-mpi \
	octopus-1d octopus-1d-mpi \
	octopus-2d octopus-2d-mpi\
	octopus-pol octopus-pol-mpi \
	zoctopus zoctopus-mpi \
	zoctopus-1d zoctopus-1d-mpi \
	zoctopus-2d zoctopus-2d-mpi \
	zoctopus-pol zoctopus-pol-mpi

octopus_SOURCES = liboct.F90 global.F90 io.F90 pulpo.F90 units.F90 math.F90 \
	spline.F90 fft.F90 \
	vxc.F90 atomic.F90 logrid.F90 tm.F90 hgh.F90 ps.F90 specie.F90 \
	atom.F90 mesh.F90 hartree.F90 out.F90 states.F90 xc.F90 \
	systm.F90 lasers.F90 h.F90 \
	eigen.F90 mix.F90 lcao.F90 scf.F90 unocc.F90 static_pol.F90 \
	geom_opt.F90 phonons.F90 pes.F90 td.F90 opt_control.F90 run.F90 main.F90

noinst_HEADERS = \
	real.F90 complex.F90 undef.F90 mesh.F90 mesh_inc.F90 mesh_create.F90\
	hartree1D.F90 hartree2D.F90 hartree3D.F90 ps3D.F90 \
	specie1D.F90 specie3D.F90\
	xc_pot.F90 xc_LDA.F90 xc_GGA.F90 xc_KLI.F90 xc_KLI_x.F90 \
	xc_KLI_SIC.F90 xc_HJU.F90 \
	h_inc.F90 h_external_pot.F90 h_forces.F90 h_so.F90 \
	eigen_cg1.F90 eigen_cg2.F90 \
	states_inc.F90 states_kpoints.F90\
	pes_rc.F90 pes_mask.F90 \
	td_init.F90 td_rti.F90 td_exp.F90 td_write.F90 td_calc.F90\
	mix_inc.F90 out_inc.F90 dx_inc.F90

octopus_LDADD = ../liboct/liboct.a @LIB_FFTW@ @LIB_LAPACK@ @LIB_BLAS@ @GSL_LIBS@

# this has to be repeated...
octopus_mpi_SOURCES     = $(octopus_SOURCES)
octopus_mpi_LDADD       = $(octopus_LDADD) @LIB_MPI@
octopus_1d_SOURCES      = $(octopus_SOURCES)
octopus_1d_LDADD        = $(octopus_LDADD)
octopus_1d_mpi_SOURCES  = $(octopus_SOURCES)
octopus_1d_mpi_LDADD    = $(octopus_mpi_LDADD)
octopus_2d_SOURCES      = $(octopus_SOURCES)
octopus_2d_LDADD        = $(octopus_LDADD)
octopus_2d_mpi_SOURCES  = $(octopus_SOURCES)
octopus_2d_mpi_LDADD    = $(octopus_mpi_LDADD)
octopus_pol_SOURCES     = $(octopus_SOURCES)
octopus_pol_LDADD       = $(octopus_LDADD)
octopus_pol_mpi_SOURCES = $(octopus_SOURCES)
octopus_pol_mpi_LDADD   = $(octopus_mpi_LDADD)
zoctopus_SOURCES        = $(octopus_SOURCES)
zoctopus_LDADD          = $(octopus_LDADD)
zoctopus_mpi_SOURCES     = $(octopus_SOURCES)
zoctopus_mpi_LDADD       = $(octopus_mpi_LDADD)
zoctopus_1d_SOURCES      = $(octopus_SOURCES)
zoctopus_1d_LDADD        = $(octopus_LDADD)
zoctopus_1d_mpi_SOURCES  = $(octopus_SOURCES)
zoctopus_1d_mpi_LDADD    = $(octopus_mpi_LDADD)
zoctopus_2d_SOURCES      = $(octopus_SOURCES)
zoctopus_2d_LDADD        = $(octopus_LDADD)
zoctopus_2d_mpi_SOURCES  = $(octopus_SOURCES)
zoctopus_2d_mpi_LDADD    = $(octopus_mpi_LDADD)
zoctopus_pol_SOURCES     = $(octopus_SOURCES)
zoctopus_pol_LDADD       = $(octopus_LDADD)
zoctopus_pol_mpi_SOURCES = $(octopus_SOURCES)
zoctopus_pol_mpi_LDADD   = $(octopus_mpi_LDADD)


# now the utility programs
strength_function_SOURCES = global.F90 io.F90 liboct.F90 units.F90 spectrum.F90 strength_function.F90
strength_function_LDADD   = ../liboct/liboct.a @GSL_LIBS@
hs_from_mult_SOURCES      = global.F90 io.F90 liboct.F90 units.F90 spectrum.F90 hs_from_mult.F90
hs_from_mult_LDADD        = ../liboct/liboct.a @GSL_LIBS@
hs_from_acc_SOURCES       = global.F90 io.F90 liboct.F90 units.F90 spectrum.F90 hs_from_acc.F90
hs_from_acc_LDADD         = ../liboct/liboct.a @GSL_LIBS@
function_to_dx_SOURCES = global.F90 io.F90 liboct.F90 units.F90 out.F90 spectrum.F90 specie.F90 atom.F90 \
                         ps.F90 spline.F90 math.F90 fft.F90 mesh.F90 states.F90 systm.F90 logrid.F90\
                         tm.F90 hgh.F90 atomic.F90 vxc.F90 mix.F90 dx.F90 function_to_dx.F90
function_to_dx_LDADD   = ../liboct/liboct.a @LIB_FFTW@ @LIB_LAPACK@ @LIB_BLAS@ @GSL_LIBS@
nbo2xyz_SOURCES = global.F90 io.F90 liboct.F90 units.F90 atom.F90 nbo2xyz.F90
nbo2xyz_LDADD   = ../liboct/liboct.a @GSL_LIBS@
phonon_spectrum_SOURCES = global.F90 io.F90 liboct.F90 units.F90 phonon_spectrum.F90
phonon_spectrum_LDADD = ../liboct/liboct.a @GSL_LIBS@
excitations_SOURCES = liboct.F90 global.F90 io.F90 units.F90 math.F90 spline.F90 \
                      fft.F90 vxc.F90 atomic.F90 logrid.F90 tm.F90 hgh.F90 ps.F90 specie.F90 \
                      atom.F90 mesh.F90 hartree.F90 out.F90 states.F90 systm.F90 mix.F90 \
                      linear.F90 excitations.F90
excitations_LDADD = ../liboct/liboct.a @LIB_FFTW@ @LIB_LAPACK@ @LIB_BLAS@ @GSL_LIBS@
broad_SOURCES = liboct.F90 global.F90 io.F90 units.F90 broad.F90
broad_LDADD = ../liboct/liboct.a @GSL_LIBS@

# cleaning
CLEANFILES = *~ *.bak *.mod *.il *.d *.pc* config_F90.h
DISTCLEANFILES = Makefile.in dirs.h

# how to compile a F90 file :)
.F90.o:
	@F90CPP@ -DSHARE_OCTOPUS='"${pkgdatadir}"' $< > $*.f90
#	sed 's/^#/!/' $*.f90 >lixo; mv -f lixo $*.f90
	@F90@ @F90FLAGS@ -c -o $@ $*.f90
	@rm -f $*.f90

# how to generate config_F90.h
liboct.F90: config_F90.h
global.F90: config_F90.h
config_F90.h: config.h
	sed 's%/\*%!%g; s%\*/%%g' config.h >config_F90.h
../liboct/liboct.a:
	(cd ../liboct; $(MAKE))

# dependencies. Please update consistently
mesh.o: mesh_inc.F90 mesh_create.F90
specie.o: specie1D.F90 specie3D.F90
ps.o: ps3D.F90
hartree.o: hartree1D.F90 hartree2D.F90 hartree3D.F90
xc.o: xc_pot.F90 xc_LDA.F90 xc_GGA.F90 xc_KLI.F90 xc_KLI_x.F90
xc.o: xc_KLI_SIC.F90 xc_HJU.F90
h.o: h_inc.F90 h_external_pot.F90 h_forces.F90 h_so.F90
eigen.o: eigen_cg1.F90 eigen_cg2.F90
states.o: states_inc.F90
pes.o: pes_rc.F90 pes_mask.F90
td.o: td_init.F90 td_rti.F90 td_exp.F90 td_write.F90 td_calc.F90
mix.o: mix_inc.F90
out.o: out_inc.F90
dx.o: dx_inc.F90
