EXTRA_PROGRAMS = oct-sf oct-rsf oct-hs-mult oct-hs-acc \
	oct-broad oct-ph oct-xyz-anim \
	oct-make-st oct-center-geom

bin_PROGRAMS = octopus oct-excite oct-atoms-magnet @UTILITY_PROGRAMS@

check_PROGRAMS = oct-test

core_SRCS = liboct_parser.F90 liboct.F90 io.F90 string.F90 global.F90 units.F90 \
	blas.F90 linalg_basic.F90 linalg_adv.F90 math.F90 spline.F90 fft.F90 \
	libxc.F90 logrid.F90 atomic.F90 tm.F90 hgh.F90 ps.F90 \
	specie.F90 xyz_file.F90 geometry.F90 \
	curv_gygi.F90 curv_briggs.F90 curv_modine.F90 curvlinear.F90 \
	mesh.F90 nl_operator.F90 \
	stencil_star.F90 stencil_variational.F90 stencil_cube.F90 derivatives.F90 \
	cf.F90 mf.F90 f.F90 \
	poisson_corrections.F90 poisson_cg.F90 poisson.F90 \
	out.F90 crystal.F90 states.F90 xc_functl.F90 xc.F90 mix.F90 restart.F90

octopus_SOURCES = $(core_SRCS) \
	lasers.F90 epot.F90 h.F90 systm.F90 \
	eigen.F90 lcao.F90 pes.F90 \
	scf.F90 gs.F90 unocc.F90 \
	static_pol.F90 geom_opt.F90 phonons.F90 \
	td_exp_split.F90 td_exp.F90 td_rti.F90 td.F90 opt_control.F90 \
	linear_response.F90 static_pol_lr.F90 pulpo.F90 \
	run.F90 main.F90

noinst_HEADERS = \
	global.h \
	real.F90 complex.F90 undef.F90 \
	fftw2.F90 fftw3.F90 \
	mesh_create.F90 derivatives_inc.F90 cf_inc.F90 mf_inc.F90 f_inc.F90 \
	poisson1D.F90 poisson2D.F90 poisson3D.F90 \
	xc_pot.F90 xc_fxc.F90 xc_KLI.F90 xc_OEP_x.F90 xc_OEP_SIC.F90 \
	epot_inc.F90 \
	h_inc.F90 h_xc_OEP.F90 h_so.F90 \
	eigen_cg2.F90 eigen_cg3.F90 eigen_plan.F90 \
	crystal.F90 restart_inc.F90 states_inc.F90 states_kpoints.F90 \
	linear_reponse_inc.F90 \
	pes_rc.F90 pes_mask.F90 \
	td_init.F90 td_write.F90 td_calc.F90\
	out_inc.F90 math_cg_inc.F90 expokit_inc.F90 math_inc.F90 \
	linalg_basic_blas.F90 linalg_basic_int.F90 linalg_adv_lapack.F90

core_LIBS = ../liboct/liboct.a ../liboct_parser/liboct_parser.a ../libxc/libxc.a \
  @GSL_LIBS@ @FCEXTRALIBS@
octopus_LDADD = @LIBS_FFT@ @LIBS_TRLAN@ @LIBS_LAPACK@ @LIBS_BLAS@ @LIBS_NETCDF@ \
	@LIBS_MPI@ $(core_LIBS) @LIBS_METIS@

# now the utility programs
oct_sf_SOURCES = global.F90 io.F90 liboct_parser.F90 liboct.F90 units.F90 spectrum.F90 strength_function.F90
oct_sf_LDADD   = $(core_LIBS)

oct_rsf_SOURCES = global.F90 io.F90 liboct_parser.F90 liboct.F90 units.F90 spectrum.F90 rotational_strength.F90
oct_rsf_LDADD   = $(core_LIBS)

oct_hs_mult_SOURCES = global.F90 io.F90 liboct_parser.F90 liboct.F90 units.F90 spectrum.F90 hs_from_mult.F90
oct_hs_mult_LDADD   = $(core_LIBS)

oct_hs_acc_SOURCES  = global.F90 io.F90 liboct_parser.F90 liboct.F90 units.F90 spectrum.F90 hs_from_acc.F90
oct_hs_acc_LDADD    = $(core_LIBS)

oct_excite_SOURCES = $(core_SRCS) linear.F90 lasers.F90 epot.F90 h.F90 systm.F90 excitations.F90
oct_excite_LDADD = $(octopus_LDADD)

oct_xyz_anim_SOURCES = $(core_SRCS) xyzanim.F90
oct_xyz_anim_LDADD   = $(octopus_LDADD)

oct_broad_SOURCES = liboct_parser.F90 liboct.F90 global.F90 io.F90 units.F90 broad.F90
oct_broad_LDADD = $(core_LIBS)

oct_ph_SOURCES = liboct_parser.F90 liboct.F90 global.F90 io.F90 units.F90 phonon_spectrum.F90
oct_ph_LDADD = $(core_LIBS)

oct_make_st_SOURCES = $(core_SRCS) lasers.F90 epot.F90 h.F90 systm.F90 make_st.F90
oct_make_st_LDADD = $(octopus_LDADD)

oct_center_geom_SOURCES = $(core_SRCS) xyz_adjust.F90 centergeom.F90
oct_center_geom_LDADD = $(octopus_LDADD)

oct_atoms_magnet_SOURCES = $(core_SRCS) atoms_magnet.F90
oct_atoms_magnet_LDADD = $(octopus_LDADD)

oct_test_SOURCES = $(core_SRCS) lasers.F90 epot.F90 h.F90 systm.F90 \
	eigen.F90 lcao.F90 pes.F90 \
	scf.F90 gs.F90 unocc.F90 \
	static_pol.F90 geom_opt.F90 phonons.F90 \
	td_exp_split.F90 td_exp.F90 td_rti.F90 td.F90 opt_control.F90 \
	linear_response.F90 static_pol_lr.F90 pulpo.F90 \
	run.F90 oct_test.F90
oct_test_LDADD = $(octopus_LDADD)

# cleaning
CLEANFILES = *~ *.bak *.mod *.il *.d *.pc* ifc* config_F90.h
DISTCLEANFILES = Makefile.in dirs.h

# ctags
CTAGS=ctags-exuberant

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
global.o: global.F90
	@FCCPP@ @CPPFLAGS@ -DSHARE_OCTOPUS='"${pkgdatadir}"' \
		-DOCTOPUS_VERSION='"${VERSION}"' -DBUILD_TIME='"${DATE}"' $< > $*_oct.f90
	@FC@ @FCFLAGS@ -c @FCFLAGS_f90@ -o $@ $*_oct.f90
	@rm -f $*_oct.f90

# how to generate config_F90.h
liboct_parser.F90: config_F90.h
global.F90: config_F90.h
config_F90.h: config.h
	sed 's%/\*%!%g; s%\*/%%g' config.h >config_F90.h
../liboct/liboct.a:
	(cd ../liboct; $(MAKE))


# This is to generate the dependencies. It makes use of the makedepf90 program of 
# Erik Edelmann (http://www.helsinki.fi/~eedelman/makedepf90.html). It should be run
# from time to time (just typing make depend within the src directory), substituting
# the dependencies below by the resulting new dependencies. (For normal compilation
# of octopus, there is absolutely no need of installing makedepf90).
depend:
	@makedepf90 $(core_SRCS) $(octopus_SOURCES) $(oct_sf_SOURCES) $(oct_rsf_SOURCES) $(oct_hs_mult_SOURCES) $(oct_hs_acc_SOURCES) $(oct_excite_SOURCES) $(oct_xyz_anim_SOURCES) $(oct_broad_SOURCES) $(oct_ph_SOURCES) $(oct_center_geom_SOURCES) $(oct_atoms_magnet_SOURCES)


# Dependencies. 
linear_response.o: linear_response_inc.F90
#liboct_parser.o : liboct_parser.F90 complex.F90 global.h 
#io.o : io.F90 liboct_parser.o 
#liboct.o : liboct.F90 complex.F90 global.h 
#global.o : global.F90 complex.F90 global.h io.o liboct.o liboct_parser.o 
#units.o : units.F90 complex.F90 global.h liboct_parser.o global.o 
#blas.o : blas.F90 complex.F90 global.h 
#linalg_basic.o : linalg_basic.F90 linalg_basic_blas.F90 complex.F90 global.h blas.o 
#linalg_adv.o : linalg_adv.F90 complex.F90 global.h global.o 
#math.o : math.F90 real.F90 math_cg_inc.F90 math_inc.F90 undef.F90 complex.F90 global.h linalg_adv.o linalg_basic.o global.o 
#spline.o : spline.F90 complex.F90 global.h 
#fft.o : fft.F90 fftw3.F90 fftw2.F90 complex.F90 global.h linalg_basic.o global.o 
#vxc.o : vxc.F90 complex.F90 global.h global.o 
#logrid.o : logrid.F90 complex.F90 global.h global.o 
#atomic.o : atomic.F90 complex.F90 global.h logrid.o vxc.o global.o 
#tm.o : tm.F90 complex.F90 global.h logrid.o io.o atomic.o global.o 
#hgh.o : hgh.F90 complex.F90 global.h logrid.o atomic.o io.o global.o 
#ps.o : ps.F90 complex.F90 global.h hgh.o tm.o atomic.o logrid.o spline.o io.o liboct.o global.o 
#specie.o : specie.F90 complex.F90 global.h math.o ps.o units.o io.o spline.o liboct_parser.o global.o 
#xyz_file.o : xyz_file.F90 complex.F90 global.h units.o io.o liboct_parser.o 
#geometry.o : geometry.F90 complex.F90 global.h xyz_file.o specie.o io.o liboct.o liboct_parser.o global.o 
#curv_gygi.o : curv_gygi.F90 complex.F90 global.h linalg_adv.o geometry.o liboct_parser.o 
#curv_briggs.o : curv_briggs.F90 complex.F90 global.h linalg_adv.o geometry.o liboct_parser.o 
#curv_modine.o : curv_modine.F90 complex.F90 global.h linalg_adv.o geometry.o liboct_parser.o 
#curvlinear.o : curvlinear.F90 complex.F90 global.h curv_modine.o curv_briggs.o curv_gygi.o geometry.o liboct_parser.o 
mesh.o : mesh.F90 mesh_create.F90 complex.F90 global.h linalg_adv.o curvlinear.o geometry.o math.o units.o liboct_parser.o global.o 
#nl_operator.o : nl_operator.F90 complex.F90 global.h mesh.o 
#stencil_star.o : stencil_star.F90 complex.F90 global.h nl_operator.o math.o global.o 
#stencil_variational.o : stencil_variational.F90 complex.F90 global.h nl_operator.o global.o 
#stencil_cube.o : stencil_cube.F90 complex.F90 global.h global.o 
#derivatives.o : derivatives.F90 derivatives_inc.F90 real.F90 undef.F90 complex.F90 global.h stencil_cube.o stencil_variational.o stencil_star.o nl_operator.o mesh.o global.o 
#cf.o : cf.F90 cf_inc.F90 real.F90 undef.F90 complex.F90 global.h fft.o mesh.o global.o 
#mf.o : mf.F90 mf_inc.F90 real.F90 undef.F90 complex.F90 global.h cf.o mesh.o linalg_basic.o global.o 
#f.o : f.F90 f_inc.F90 real.F90 undef.F90 complex.F90 global.h fft.o derivatives.o cf.o mf.o linalg_basic.o liboct_parser.o global.o 
#poisson.o : poisson.F90 poisson3D.F90 poisson2D.F90 poisson1D.F90 complex.F90 global.h math.o f.o mesh.o geometry.o liboct_parser.o global.o 
#out.o : out.F90 out_inc.F90 real.F90 undef.F90 complex.F90 global.h f.o mesh.o units.o io.o liboct_parser.o liboct.o global.o 
#crystal.o : crystal.F90 complex.F90 global.h units.o io.o global.o 
#states.o : states.F90 states_inc.F90 real.F90 undef.F90 states_kpoints.F90 complex.F90 global.h crystal.o geometry.o out.o mf.o f.o mesh.o math.o linalg_adv.o linalg_basic.o io.o liboct_parser.o global.o 
xc.o : xc.F90 xc_OEP_SIC.F90 xc_OEP_x.F90 xc_KLI.F90 xc_pot.F90
#mix.o : mix.F90 complex.F90 global.h linalg_adv.o linalg_basic.o liboct_parser.o global.o 
#restart.o : restart.F90 complex.F90 global.h out.o mesh.o states.o io.o global.o 
#lasers.o : lasers.F90 complex.F90 global.h mesh.o spline.o liboct_parser.o io.o global.o 
#epot.o : epot.F90 epot_inc.F90 real.F90 undef.F90 complex.F90 global.h lasers.o states.o geometry.o specie.o ps.o logrid.o mesh.o liboct_parser.o global.o 
#systm.o : systm.F90 complex.F90 global.h states.o geometry.o mesh.o global.o 
#h.o : h.F90 h_so.F90 h_xc_OEP.F90 h_inc.F90 real.F90 undef.F90 complex.F90 global.h out.o xc.o poisson.o epot.o states.o geometry.o mf.o mesh.o linalg_basic.o liboct_parser.o global.o 
#eigen.o : eigen.F90 eigen_plan.F90 eigen_cg2.F90 complex.F90 global.h h.o states.o mesh.o stencil_star.o nl_operator.o io.o linalg_adv.o linalg_basic.o liboct_parser.o global.o 
#lcao.o : lcao.F90 complex.F90 global.h h.o systm.o geometry.o specie.o mesh.o linalg_adv.o linalg_basic.o spline.o liboct.o global.o 
#pes.o : pes.F90 pes_mask.F90 pes_rc.F90 complex.F90 global.h states.o mesh.o fft.o units.o io.o liboct_parser.o global.o 
#scf.o : scf.F90 complex.F90 global.h lcao.o mix.o eigen.o h.o restart.o states.o geometry.o units.o liboct_parser.o global.o 
#gs.o : gs.F90 complex.F90 global.h scf.o lcao.o h.o systm.o global.o 
#unocc.o : unocc.F90 complex.F90 global.h eigen.o h.o restart.o systm.o states.o mesh.o liboct.o liboct_parser.o global.o 
#static_pol.o : static_pol.F90 complex.F90 global.h scf.o geometry.o states.o h.o io.o liboct.o liboct_parser.o global.o 
#geom_opt.o : geom_opt.F90 complex.F90 global.h scf.o states.o geometry.o h.o mesh.o global.o 
#phonons.o : phonons.F90 complex.F90 global.h scf.o states.o h.o geometry.o epot.o linalg_adv.o io.o liboct_parser.o liboct.o global.o 
#td_exp_split.o : td_exp_split.F90 complex.F90 global.h h.o states.o mesh.o linalg_basic.o global.o 
#td_exp.o : td_exp.F90 complex.F90 global.h td_exp_split.o mesh.o math.o linalg_basic.o blas.o liboct_parser.o global.o 
#td_rti.o : td_rti.F90 complex.F90 global.h td_exp.o h.o states.o mesh.o global.o 
#td.o : td.F90 td_init.F90 td_calc.F90 td_write.F90 complex.F90 global.h pes.o td_rti.o epot.o h.o restart.o states.o mesh.o liboct_parser.o liboct.o global.o 
#opt_control.o : opt_control.F90 complex.F90 global.h td_rti.o td.o h.o restart.o states.o io.o liboct_parser.o liboct.o global.o 
#pulpo.o : pulpo.F90 liboct.o global.o io.o 
#run.o : run.F90 complex.F90 global.h pulpo.o gs.o opt_control.o phonons.o geom_opt.o static_pol.o unocc.o scf.o systm.o h.o epot.o restart.o states.o lcao.o global.o 
#main.o : main.F90 complex.F90 global.h run.o liboct.o global.o 
#spectrum.o : spectrum.F90 complex.F90 global.h units.o io.o global.o 
#strength_function.o : strength_function.F90 complex.F90 global.h spectrum.o units.o liboct_parser.o global.o 
#rotational_strength.o : rotational_strength.F90 complex.F90 global.h spectrum.o units.o liboct_parser.o global.o 
#hs_from_mult.o : hs_from_mult.F90 complex.F90 global.h spectrum.o liboct_parser.o units.o global.o 
#hs_from_acc.o : hs_from_acc.F90 complex.F90 global.h spectrum.o liboct_parser.o units.o global.o 
#linear.o : linear.F90 complex.F90 global.h states.o poisson.o mf.o mesh.o linalg_adv.o io.o liboct.o global.o 
#excitations.o : excitations.F90 complex.F90 global.h linear.o systm.o restart.o states.o liboct.o liboct_parser.o global.o 
#xyzanim.o : xyzanim.F90 complex.F90 global.h geometry.o units.o io.o liboct_parser.o liboct.o global.o 
#broad.o : broad.F90 complex.F90 global.h liboct_parser.o units.o io.o global.o 
#phonon_spectrum.o : phonon_spectrum.F90 complex.F90 global.h units.o io.o global.o 
#xyz_adjust.o : xyz_adjust.F90 complex.F90 global.h geometry.o 
#centergeom.o : centergeom.F90 complex.F90 global.h xyz_adjust.o geometry.o units.o global.o 
#atoms_magnet.o : atoms_magnet.F90 complex.F90 global.h restart.o states.o out.o geometry.o mf.o mesh.o units.o liboct.o liboct_parser.o global.o 
