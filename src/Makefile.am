## Process this file with automake to produce Makefile.in

## Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
## 02111-1307, USA.
##
## $Id$

SUBDIRS = include basic math species grid

EXTRA_PROGRAMS = oct-cross-section oct-rsf oct-harmonic-spectrum \
                 oct-broad oct-ph oct-xyz-anim oct-center-geom oct-help \
                 oct-oscillator-strength

bin_PROGRAMS = octopus oct-test @UTILITY_PROGRAMS@

core_SRCS = \
	poisson_isf.F90 poisson_fft.F90 \
	poisson_corrections.F90 poisson_cg.F90 poisson_multigrid.F90 poisson.F90 \
	crystal.F90 \
	states.F90 magnetic.F90 elf.F90 states_output.F90 excited_states.F90 \
	xc_functl.F90 xc.F90 scf_tol.F90 mix.F90 spectrum.F90 \
	specie_pot.F90 tdfunction.F90

octopus_SOURCES = $(core_SRCS) \
	lasers.F90 \
	rkb_projector.F90 kb_projector.F90 hgh_projector.F90 projector.F90 \
	epot.F90 \
	h.F90 preconditioners.F90 linear_solver.F90 linear_response.F90 restart.F90 \
	xc_OEP.F90 v_ks.F90 systm.F90  xc_OEP_kernel.F90\
	eigen.F90 lcao.F90 pes.F90 \
	scf.F90 gs.F90 unocc.F90 \
	static_pol.F90 geom_opt.F90 phonons.F90 phonons_fd.F90 \
	td_exp_split.F90 td_exp.F90 td_rti.F90 td_write.F90 td.F90 filters.F90 opt_control_parameters.F90 \
	opt_control.F90 casida.F90 td_transport.F90 resp_pert.F90 sternheimer.F90 em_resp_calc.F90  \
	em_resp.F90 vdW.F90 phonons_lr.F90 pulpo.F90 run.F90 main.F90

noinst_HEADERS = \
	poisson1D.F90 poisson2D.F90 poisson3D.F90 mix_inc.F90 \
	xc_vxc.F90 xc_axc.F90 xc_fxc.F90 xc_kxc.F90 xc_KLI_inc.F90 xc_OEP_x_inc.F90 xc_OEP_SIC_inc.F90 \
	kb_projector_inc.F90 hgh_projector_inc.F90 projector_inc.F90 epot_inc.F90 \
	xc_OEP_inc.F90 h_inc.F90 \
	eigen_trlan.F90 preconditioners_inc.F90 \
	eigen_inc.F90 eigen_cg_inc.F90 eigen_plan_inc.F90 \
	eigen_arpack_inc.F90 eigen_evolution_inc.F90 \
	lcao_inc.F90 crystal.F90 restart_inc.F90 states_inc.F90 states_kpoints.F90 \
	pes_rc.F90 pes_mask.F90 \
	td_init.F90 td_calc.F90\
	linear_response_inc.F90 linear_response_out.F90 linear_solver_inc.F90\
	sternheimer_inc.F90 em_resp_inc.F90 resp_pert_inc.F90 unocc_inc.F90 casida_inc.F90 \
	excited_states_inc.F90 xc_OEP_kernel_inc.F90 em_resp_calc_inc.F90 \
	opt_control_read.F90 opt_control_aux.F90 opt_control_defstates.F90 \
	opt_control_propagation.F90 opt_control_finalcheck.F90 opt_control_iter.F90 \
	opt_control_output.F90 opt_control_tdtarget.F90

octopus_LIBS = \
	grid/libgrid.a species/libspecies.a math/libmath.a basic/libbasic.a

core_LIBS = $(octopus_LIBS) @LIBS_LAPACK@ @LIBS_BLAS@ \
  ../liboct/liboct.a \
  ../liboct_parser/liboct_parser.a \
  ../libxc/src/libxc.a \
  -L../libstring_f -lstring_f \
  @GSL_LIBS@ @GD_LIBS@ @FCEXTRALIBS@

external_LIBS = ../external_libs/expokit/libexpokit.a \
  ../external_libs/qshep/libqshep.a \
  ../external_libs/poisson_isf/libpoisson_isf.a

if COMPILE_METIS
  external_LIBS += ../external_libs/metis-4.0/libmetis.a
endif

if COMPILE_LIBNBC
  external_LIBS += ../external_libs/libnbc/libnbc.a
endif

octopus_LDADD = $(core_LIBS) @LIBS_FFT@ @LIBS_TRLAN@ @LIBS_ARPACK@ @LIBS_SPARSKIT@ \
  @LIBS_NETCDF@ $(external_LIBS) @LIBS_MPI@ @LIBS_LAPACK@ @LIBS_BLAS@

# now the utility programs
oct_cross_section_SOURCES = $(basic_SRCS) spectrum.F90 cross_section.F90
oct_cross_section_LDADD   = @LIBS_LAPACK@ @LIBS_BLAS@ $(core_LIBS) $(external_LIBS)

oct_oscillator_strength_SOURCES = getopt_f.c command_line.F90 $(basic_SRCS) spectrum.F90 oscillator_strength.F90
oct_oscillator_strength_LDADD   = @LIBS_LAPACK@ @LIBS_BLAS@ $(core_LIBS) $(external_LIBS)

oct_rsf_SOURCES = $(basic_SRCS) spectrum.F90 rotational_strength.F90
oct_rsf_LDADD   = @LIBS_LAPACK@ @LIBS_BLAS@ $(core_LIBS) $(external_LIBS)

oct_harmonic_spectrum_SOURCES = $(basic_SRCS) spectrum.F90 hs.F90
oct_harmonic_spectrum_LDADD   = @LIBS_LAPACK@ @LIBS_BLAS@ $(core_LIBS) $(external_LIBS)

oct_xyz_anim_SOURCES = $(core_SRCS) xyzanim.F90
oct_xyz_anim_LDADD   = $(octopus_LDADD)

oct_broad_SOURCES = $(basic_SRCS) broad.F90
oct_broad_LDADD =  $(core_LIBS) $(external_LIBS)

oct_ph_SOURCES = $(basic_SRCS) phonon_spectrum.F90
oct_ph_LDADD =  $(core_LIBS) $(external_LIBS)

oct_center_geom_SOURCES = $(core_SRCS) centergeom.F90
oct_center_geom_LDADD = $(octopus_LDADD)

oct_test_SOURCES = $(core_SRCS) lasers.F90 rkb_projector.F90 kb_projector.F90 \
	hgh_projector.F90 projector.F90 epot.F90 h.F90 preconditioners.F90 linear_solver.F90 \
	linear_response.F90 restart.F90 xc_OEP.F90 v_ks.F90 \
	systm.F90 test.F90
oct_test_LDADD = $(octopus_LDADD)

oct_help_SOURCES = getopt_f.c command_line.F90 $(basic_SRCS) help.F90
oct_help_LDADD =  $(core_LIBS)

AM_FCFLAGS = \
	@F90_MODULE_FLAG@$(srcdir)/basic   \
	@F90_MODULE_FLAG@$(srcdir)/math    \
	@F90_MODULE_FLAG@$(srcdir)/species \
	@F90_MODULE_FLAG@$(srcdir)/grid \
	@F90_MODULE_FLAG@$(top_builddir)/libxc/src/ \
	@F90_MODULE_FLAG@$(top_builddir)/external_libs/qshep

AM_CPPFLAGS = -I$(srcdir)/include -I$(top_srcdir)/libstring_f

include common-rules.make

# if it exists include user defined Makefile.local (which is not included in the octopus package)
-include Makefile.local

# cleaning
DISTCLEANFILES = dirs.h


../liboct/liboct.a:
	(cd ../liboct; $(MAKE))
basic/libbasic.a:
	(cd basic; $(MAKE))
math/libmath.a:
	(cd math; $(MAKE))
species/libspecies.a:
	(cd species; $(MAKE))


# This is to generate the dependencies. It makes use of the makedepf90 program of 
# Erik Edelmann (http://www.helsinki.fi/~eedelman/makedepf90.html). It should be run
# from time to time (just typing make depend within the src directory), substituting
# the dependencies below by the resulting new dependencies. (For normal compilation
# of octopus, there is absolutely no need of installing makedepf90).
depend:
	@makedepf90 $(core_SRCS) $(octopus_SOURCES) $(oct_cross_section_SOURCES) $(oct_rsf_SOURCES) $(oct_harmonic_spectrum_SOURCES) $(oct_xyz_anim_SOURCES) $(oct_broad_SOURCES) $(oct_ph_SOURCES) $(oct_center_geom_SOURCES) $(oct_atoms_magnet_SOURCES) $(oct_help_SOURCES)


# Dependencies (only includes taken into account)
casida.o:  casida_inc.F90
eigen.o:  eigen_inc.F90 eigen_cg_inc.F90 eigen_trlan.F90 eigen_plan_inc.F90 eigen_arpack_inc.F90 eigen_evolution_inc.F90
sternheimer.o:  sternheimer_inc.F90
em_resp_calc.o:  em_resp_calc_inc.F90
em_resp.o :  em_resp_inc.F90
epot.o:  epot_inc.F90
excited_states.o:  excited_states_inc.F90
fft.o:  fftw2.F90 fftw3.F90
h.o:  h_inc.F90
hgh_projector.o:  hgh_projector_inc.F90
kb_projector.o:  kb_projector_inc.F90
lcao.o:  lcao_inc.F90
linear_solver.o:  linear_solver_inc.F90
linear_response.o:  linear_response_inc.F90 linear_response_out.F90
mix.o:  mix_inc.F90
opt_control.o:  opt_control_read.F90 opt_control_aux.F90 opt_control_defstates.F90 \
               opt_control_propagation.F90 opt_control_finalcheck.F90 opt_control_iter.F90 \
               opt_control_output.F90 opt_control_tdtarget.F90
pes.o:  pes_rc.F90 pes_mask.F90
poisson.o:  poisson1D.F90 poisson2D.F90 poisson3D.F90 poisson_multigrid.F90
projector.o:  projector_inc.F90
restart.o:  restart_inc.F90
preconditioners.o:  preconditioners_inc.F90
states.o:  states_kpoints.F90 states_inc.F90
em_resp.o:  em_resp_inc.F90
resp_pert.o: resp_pert_inc.F90
td.o:  td_init.F90
td_write.o: td_calc.F90
unocc.o:  unocc_inc.F90
xc.o:  xc_vxc.F90 xc_axc.F90 xc_fxc.F90 xc_kxc.F90
xc_OEP.o:  xc_KLI_inc.F90 xc_OEP_x_inc.F90 xc_OEP_SIC_inc.F90 xc_OEP_inc.F90
xc_OEP_kernel.o:  xc_OEP_kernel_inc.F90