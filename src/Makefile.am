EXTRA_PROGRAMS = oct-sf oct-rsf oct-hs-mult oct-hs-acc \
	oct-excite oct-broad oct-ph oct-xyz-anim \
	oct-choose-st oct-make-st oct-center-geom \
	oct-atoms-magnet

bin_PROGRAMS = octopus @UTILITY_PROGRAMS@

core_SRCS = liboct_parser.F90 liboct.F90 global.F90 io.F90 units.F90 \
	linalg_basic.F90 linalg_adv.F90 math.F90 spline.F90 fft.F90 \
	vxc.F90 logrid.F90 atomic.F90 tm.F90 hgh.F90 ps.F90 \
	specie.F90 xyz_file.F90 geometry.F90 \
	mesh.F90 cf.F90 mf.F90 f.F90 \
	poisson.F90 out.F90 crystal.F90 states.F90 xc.F90 mix.F90 restart.F90

octopus_SOURCES = $(core_SRCS) \
	lasers.F90 epot.F90 systm.F90 h.F90 \
	eigen.F90 lcao.F90 pes.F90 \
	scf.F90 unocc.F90 \
	static_pol.F90 geom_opt.F90 phonons.F90 \
	td_exp_split.F90 td_exp.F90 td_rti.F90 td.F90 \
	opt_control.F90 pulpo.F90 \
	run.F90 main.F90

noinst_HEADERS = \
	global.h \
	real.F90 complex.F90 undef.F90 \
	fftw2.F90 fftw3.F90 \
	mesh_create.F90 cf_inc.F90 mf_inc.F90 f_inc.F90 \
	poisson1D.F90 poisson2D.F90 poisson3D.F90 \
	specie1D.F90 specie3D.F90\
	xc_pot.F90 xc_LDA.F90 xc_GGA.F90 xc_KLI.F90 xc_OEP_x.F90 xc_OEP_SIC.F90 \
	epot_inc.F90 \
	h_inc.F90 h_xc_OEP.F90 h_so.F90 \
	eigen_cg2.F90 eigen_cg3.F90 eigen_plan.F90\
	crystal.F90 states_inc.F90 states_kpoints.F90\
	pes_rc.F90 pes_mask.F90 \
	td_init.F90 td_write.F90 td_calc.F90\
	out_inc.F90 math_inc.F90\
	linalg_basic_blas.F90 linalg_basic_int.F90 linalg_adv_lapack.F90

core_LIBS = ../liboct/liboct.a ../liboct/liboct_parser.a @GSL_LIBS@ @LDF90FLAGS@
octopus_LDADD = @LIBS_FFTW@ @LIBS_TRLAN@ @LIBS_LAPACK@ @LIBS_BLAS@ @LIBS_NETCDF@ \
	@LIBS_MPI@ $(core_LIBS)

# now the utility programs
oct_sf_SOURCES = global.F90 io.F90 liboct_parser.F90 liboct.F90 units.F90 spectrum.F90 strength_function.F90
oct_sf_LDADD   = $(core_LIBS)

oct_rsf_SOURCES = global.F90 io.F90 liboct_parser.F90 liboct.F90 units.F90 spectrum.F90 rotational_strength.F90
oct_rsf_LDADD   = $(core_LIBS)

oct_hs_mult_SOURCES = global.F90 io.F90 liboct_parser.F90 liboct.F90 units.F90 spectrum.F90 hs_from_mult.F90
oct_hs_mult_LDADD   = $(core_LIBS)

oct_hs_acc_SOURCES  = global.F90 io.F90 liboct_parser.F90 liboct.F90 units.F90 spectrum.F90 hs_from_acc.F90
oct_hs_acc_LDADD    = $(core_LIBS)

oct_excite_SOURCES = $(core_SRCS) linear.F90 systm.F90 excitations.F90
oct_excite_LDADD = $(octopus_LDADD)

oct_xyz_anim_SOURCES = $(core_SRCS) xyzanim.F90
oct_xyz_anim_LDADD   = $(octopus_LDADD)

oct_broad_SOURCES = liboct_parser.F90 liboct.F90 global.F90 io.F90 units.F90 broad.F90
oct_broad_LDADD = $(core_LIBS)

oct_ph_SOURCES = liboct_parser.F90 liboct.F90 global.F90 io.F90 units.F90 phonon_spectrum.F90
oct_ph_LDADD = $(core_LIBS)

oct_choose_st_SOURCES = $(core_SRCS) systm.F90 choose_st.F90
oct_choose_st_LDADD = $(octopus_LDADD)

oct_make_st_SOURCES = $(core_SRCS) systm.F90 make_st.F90
oct_make_st_LDADD = $(octopus_LDADD)

oct_center_geom_SOURCES = $(core_SRCS) xyz_adjust.F90 centergeom.F90
oct_center_geom_LDADD = $(octopus_LDADD)

oct_atoms_magnet_SOURCES = $(core_SRCS) atoms_magnet.F90
oct_atoms_magnet_LDADD = $(octopus_LDADD)

# cleaning
CLEANFILES = *~ *.bak *.mod *.il *.d *.pc* ifc* config_F90.h
DISTCLEANFILES = Makefile.in dirs.h

# how to compile a F90 file :)
.F90.o:
	@FCCPP@ $< > $*.f90
	@if [ "x@DEBUG@" = "x" ]; then \
		cat $*.f90 | grep -v pop_sub | grep -v push_sub >$*.f91; \
		mv -f $*.f91 $*.f90; \
	fi
	@FC@ @FCFLAGS@ -c @FCFLAGS_f90@ -o $@ $*.f90
	@rm -f $*.f90

# how to compile main.f90 (have to define some extra variables here)
DATE := $(shell date)
global.o: global.F90
	@FCCPP@ -DSHARE_OCTOPUS='"${pkgdatadir}"' \
		-DOCTOPUS_VERSION='"${VERSION}"' -DBUILD_TIME='"${DATE}"' $< > $*.f90
	@FC@ @FCFLAGS@ -c @FCFLAGS_f90@ -o $@ $*.f90
	@rm -f $*.f90

# how to generate config_F90.h
liboct_parser.F90: config_F90.h
global.F90: config_F90.h
config_F90.h: config.h
	sed 's%/\*%!%g; s%\*/%%g' config.h >config_F90.h
../liboct/liboct.a:
	(cd ../liboct; $(MAKE))

# dependencies. Please update consistently
mesh.o: mesh_create.F90
cf.o: cf_inc.F90
mf.o: mf_inc.F90
f.o: f_inc.F90
specie.o: specie1D.F90 specie3D.F90
poisson.o: poisson1D.F90 poisson2D.F90 poisson3D.F90
xc.o: xc_pot.F90 xc_LDA.F90 xc_GGA.F90
xc.o: xc_KLI.F90 xc_OEP_x.F90 xc_OEP_SIC.F90
epot.o: epot_inc.F90
h.o: h_inc.F90 h_xc_OEP.F90 h_so.F90
eigen.o: eigen_cg2.F90 eigen_cg3.F90 eigen_plan.F90
states.o: states_inc.F90 states_kpoints.F90
pes.o: pes_rc.F90 pes_mask.F90
td.o: td_init.F90 td_write.F90 td_calc.F90
out.o: out_inc.F90
math.o: math_inc.F90
fft.o: fftw2.F90 fftw3.F90
linalg_basic.o : linalg_basic_blas.F90 linalg_basic_int.F90
linalg_adv.o : linalg_adv_lapack.F90
