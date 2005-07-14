# $Id$

EXTRA_PROGRAMS = oct-sf oct-rsf oct-hs-mult oct-hs-acc \
	oct-broad oct-ph oct-xyz-anim \
	oct-make-st oct-center-geom

bin_PROGRAMS = octopus @UTILITY_PROGRAMS@

check_PROGRAMS = oct-test

basic_SRCS = liboct_parser.F90 liboct.F90 string.F90 varinfo.F90 global.F90 \
	     messages.F90 syslabels.F90 io.F90 units.F90

core_SRCS = $(basic_SRCS) \
	blas.F90 linalg_basic.F90 linalg_adv.F90 math.F90 spline.F90 sparskit.F90 \
	libxc.F90 logrid.F90 atomic.F90 tm.F90 hgh.F90 ps.F90 \
	specie.F90 xyz_file.F90 geometry.F90 simul_box.F90 \
	curv_gygi.F90 curv_briggs.F90 curv_modine.F90 curvlinear.F90 \
	fft.F90 mesh.F90 nl_operator.F90 \
	stencil_star.F90 stencil_variational.F90 stencil_cube.F90 stencil_starplus.F90 derivatives.F90 \
	cf.F90 mf.F90 f.F90 \
	grid.F90 \
	out.F90 poisson_corrections.F90 poisson_cg.F90 poisson.F90 \
	crystal.F90 states.F90 xc_functl.F90 xc.F90 mix.F90 restart.F90

octopus_SOURCES = $(core_SRCS) \
	lasers.F90 epot.F90 h.F90 linear_response.F90 xc_OEP.F90 v_ks.F90 systm.F90 \
	eigen.F90 lcao.F90 pes.F90 \
	scf.F90 gs.F90 unocc.F90 \
	static_pol.F90 geom_opt.F90 phonons.F90 \
	td_exp_split.F90 td_exp.F90 td_rti.F90 td.F90 opt_control.F90 \
	casida.F90 wave_matching.F90 static_pol_lr.F90 pulpo.F90 \
	run.F90 main.F90

noinst_HEADERS = \
	global.h \
	real.F90 complex.F90 undef.F90 \
	fftw2.F90 fftw3.F90 \
	mesh_create.F90 derivatives_inc.F90 cf_inc.F90 mf_inc.F90 f_inc.F90 \
	poisson1D.F90 poisson2D.F90 poisson3D.F90 \
	xc_vxc.F90 xc_axc.F90 xc_fxc.F90 xc_KLI.F90 xc_OEP_x.F90 xc_OEP_SIC.F90 \
	epot_inc.F90 \
	v_ks_inc.F90 xc_OEP_inc.F90 h_inc.F90 h_so.F90 \
	eigen_cg2.F90 eigen_cg3.F90 eigen_plan.F90 eigen_arpack.F90 eigen_evolution.F90 \
	crystal.F90 restart_inc.F90 states_inc.F90 states_kpoints.F90 \
	linear_response_inc.F90 \
	pes_rc.F90 pes_mask.F90 \
	td_init.F90 td_write.F90 td_calc.F90\
	out_inc.F90 math_cg_inc.F90 math_inc.F90 \
	linalg_basic_blas.F90 linalg_basic_int.F90 linalg_adv_lapack.F90

core_LIBS = ../liboct/liboct.a ../liboct_parser/liboct_parser.a ../libxc/libxc.a ../expokit/libexpokit.a \
  @GSL_LIBS@ @FCEXTRALIBS@
octopus_LDADD = @LIBS_FFT@ @LIBS_TRLAN@ @LIBS_ARPACK@ @LIBS_SPARSKIT@ @LIBS_LAPACK@ @LIBS_BLAS@ @LIBS_NETCDF@ \
	@LIBS_MPI@ $(core_LIBS) @LIBS_METIS@

# now the utility programs
oct_sf_SOURCES = $(basic_SRCS) spectrum.F90 strength_function.F90
oct_sf_LDADD   = $(core_LIBS)

oct_rsf_SOURCES = $(basic_SRCS) spectrum.F90 rotational_strength.F90
oct_rsf_LDADD   = $(core_LIBS)

oct_hs_mult_SOURCES = $(basic_SRCS) spectrum.F90 hs_from_mult.F90
oct_hs_mult_LDADD   = $(core_LIBS)

oct_hs_acc_SOURCES  = $(basic_SRCS) spectrum.F90 hs_from_acc.F90
oct_hs_acc_LDADD    = $(core_LIBS)

oct_xyz_anim_SOURCES = $(core_SRCS) xyzanim.F90
oct_xyz_anim_LDADD   = $(octopus_LDADD)

oct_broad_SOURCES = $(basic_SRCS) broad.F90
oct_broad_LDADD = $(core_LIBS)

oct_ph_SOURCES = $(basic_SRCS) phonon_spectrum.F90
oct_ph_LDADD = $(core_LIBS)

oct_make_st_SOURCES = $(core_SRCS) lasers.F90 epot.F90 h.F90 linear_response.F90 \
	xc_OEP.F90 v_ks.F90 systm.F90 make_st.F90
oct_make_st_LDADD = $(octopus_LDADD)

oct_center_geom_SOURCES = $(core_SRCS) xyz_adjust.F90 centergeom.F90
oct_center_geom_LDADD = $(octopus_LDADD)

oct_test_SOURCES = $(core_SRCS) lasers.F90 epot.F90 h.F90 linear_response.F90 xc_OEP.F90 v_ks.F90 systm.F90 \
	eigen.F90 lcao.F90 pes.F90 \
	scf.F90 gs.F90 unocc.F90 \
	static_pol.F90 geom_opt.F90 phonons.F90 \
	td_exp_split.F90 td_exp.F90 td_rti.F90 td.F90 opt_control.F90 \
	casida.F90 wave_matching.F90 static_pol_lr.F90 pulpo.F90 \
	run.F90 oct_test.F90

oct_test_LDADD = $(octopus_LDADD)

# cleaning
CLEANFILES = *~ *.bak *.mod *.il *.d *.pc* ifc* config_F90.h
DISTCLEANFILES = Makefile.in dirs.h

# ctags
CTAGS=ctags-exuberant -e

# if it exists include user defined Makefile.local (which is not included in the octopus package)
-include Makefile.local

# how to compile a F90 file :)
.F90.o:
	@FCCPP@ @CPPFLAGS@ $< > $*_oct.f90
	@if [ "x@DEBUG@" = "x" ]; then \
		cat $*_oct.f90 | grep -v pop_sub | grep -v push_sub >$*_oct.f91; \
		mv -f $*_oct.f91 $*_oct.f90; \
	fi
	@FC@ @FCFLAGS@ -c @FCFLAGS_f90@ -o $@ $*_oct.f90
	@rm -f $*_oct.f90

# how to compile main.f90 (have to define some extra variables here)
DATE := $(shell date)
LATEST_CVS := $(shell ../build/latest_cvs_commit.sh 2>/dev/null) 

global.o: global.F90
	@FCCPP@ @CPPFLAGS@ -DSHARE_OCTOPUS='"${pkgdatadir}"' -DLATEST_CVS='"${LATEST_CVS}"' \
		-DOCTOPUS_VERSION='"${VERSION}"' -DBUILD_TIME='"${DATE}"' $< > $*_oct.f90
	@FC@ @FCFLAGS@ -c @FCFLAGS_f90@ -o $@ $*_oct.f90
	@rm -f $*_oct.f90

# how to generate config_F90.h
liboct_parser.F90: config_F90.h
global.F90: config_F90.h
config_F90.h: ../config.h
	grep '#define' ../config.h >config_F90.h
../liboct/liboct.a:
	(cd ../liboct; $(MAKE))


# This is to generate the dependencies. It makes use of the makedepf90 program of 
# Erik Edelmann (http://www.helsinki.fi/~eedelman/makedepf90.html). It should be run
# from time to time (just typing make depend within the src directory), substituting
# the dependencies below by the resulting new dependencies. (For normal compilation
# of octopus, there is absolutely no need of installing makedepf90).
depend:
	@makedepf90 $(core_SRCS) $(octopus_SOURCES) $(oct_sf_SOURCES) $(oct_rsf_SOURCES) $(oct_hs_mult_SOURCES) $(oct_hs_acc_SOURCES) $(oct_xyz_anim_SOURCES) $(oct_broad_SOURCES) $(oct_ph_SOURCES) $(oct_center_geom_SOURCES) $(oct_atoms_magnet_SOURCES)


# Dependencies (only includes taken into account)
atomic.o: global.h
blas.o: global.h
broad.o: global.h
casida.o: global.h
centergeom.o: global.h
cf.o: global.h cf_inc.F90
crystal.o: global.h
curv_briggs.o: global.h
curv_gygi.o: global.h
curvlinear.o: global.h
curv_modine.o: global.h
derivatives.o: global.h derivatives_inc.F90
eigen.o: global.h eigen_cg2.F90 eigen_cg3.F90 eigen_plan.F90 eigen_arpack.F90 eigen_evolution.F90
epot.o: global.h epot_inc.F90 epot_inc.F90
f.o: global.h f_inc.F90
fft.o: global.h fftw2.F90 fftw3.F90
geometry.o: global.h
geom_opt.o: global.h
gs.o: global.h
h.o: global.h h_inc.F90 h_so.F90
hgh.o: global.h
hs_from_acc.o: global.h
hs_from_mult.o: global.h
lasers.o: global.h
lcao.o: global.h
liboct.o: global.h
liboct_parser.o: global.h
libxc.o: global.h
linalg_adv.o: global.h linalg_adv_lapack.F90
linalg_basic.o: global.h linalg_basic_blas.F90 linalg_basic_int.F90
linear_response.o: global.h linear_response_inc.F90
logrid.o: global.h
main.o: global.h
make_st.o: global.h
math.o: global.h math_cg_inc.F90 math_inc.F90
mesh.o: global.h mesh_create.F90
mf.o: global.h mf_inc.F90
mix.o: global.h
nl_operator.o: global.h
oct_test.o: global.h
opt_control.o: global.h
out.o: global.h out_inc.F90
pes.o: global.h pes_rc.F90 pes_mask.F90
phonons.o: global.h
phonon_spectrum.o: global.h
poisson_cg.o: global.h
poisson_corrections.o: global.h
poisson.o: global.h poisson1D.F90 poisson2D.F90 poisson3D.F90
ps.o: global.h
restart.o: global.h restart_inc.F90
rotational_strength.o: global.h
run.o: global.h
scf.o: global.h
sparskit.o: global.h
specie.o: global.h
spectrum.o: global.h
spline.o: global.h
states.o: global.h states_kpoints.F90 states_inc.F90
static_pol.o: global.h
static_pol_lr.o: global.h
stencil_cube.o: global.h
stencil_star.o: global.h
stencil_variational.o: global.h
strength_function.o: global.h
systm.o: global.h
td_exp.o: global.h
td_exp_split.o: global.h
td.o: global.h td_write.F90 td_init.F90
td_rti.o: global.h
td_write.o: td_calc.F90
tm.o: global.h
units.o: global.h
unocc.o: global.h
v_ks.o: global.h v_ks_inc.F90
xc.o: global.h xc_vxc.F90 xc_axc.F90 xc_fxc.F90
xc_functl.o: global.h
xc_OEP.o: global.h xc_KLI.F90 xc_OEP_x.F90 xc_OEP_SIC.F90 xc_OEP_inc.F90
xyz_adjust.o: global.h
xyzanim.o: global.h
xyz_file.o: global.h

# manually sort out a few dependencies
varinfo.o: string.o
global.o: global.h varinfo.o string.o
messages.o: liboct.o
liboct_parser.o: global.o messages.o
syslabels: liboct_parser.o
io.o: global.o liboct_parser.o syslabels.o



