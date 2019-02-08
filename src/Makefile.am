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
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
## 02110-1301, USA.
##

include common-rules.make

noinst_HEADERS =                          \
	include/global.h                  \
	include/integer.F90               \
	include/real.F90                  \
	include/complex.F90               \
	include/real_single.F90           \
	include/complex.F90               \
	include/complex_single.F90        \
	include/undef.F90                 \
	include/vectors.h                 \
	include/fortran_types.h

include/config_F90.h: $(top_builddir)/config.h
	$(MKDIR_P) include # required for VPATH functionality
	grep '^#define' $^ > include/config_F90.h

# include/options.h and include/defaults.h rules are after the source declaration

# Fortran and C sources - keep in alphabetical order for each folder
# Sources in basic
all_f_srcs =                                   \
	basic/accel.F90                        \
	basic/blacs.F90                        \
	basic/blacs_proc_grid.F90              \
	basic/calc_mode_par.F90                \
	basic/command_line.F90                 \
	basic/comm.F90                         \
	basic/cuda.F90                         \
	basic/debug.F90                        \
	basic/distributed.F90                  \
	basic/global.F90                       \
	basic/hardware.F90                     \
	basic/ialist.F90                       \
	basic/iihash.F90                       \
	basic/io_binary_f.F90                  \
	basic/io.F90                           \
	basic/loct.F90                         \
	basic/loct_pointer.F90                 \
	basic/lookup.F90                       \
	basic/messages.F90                     \
	basic/mpi_debug.F90                    \
	basic/mpi.F90                          \
	basic/mpi_lib.F90                      \
	basic/multicomm.F90                    \
	basic/parser.F90                       \
	basic/profiling.F90                    \
	basic/sort.F90                         \
	basic/space.F90                        \
	basic/string.F90                       \
	basic/types.F90                        \
	basic/unit.F90                         \
	basic/unit_system.F90                  \
	basic/utils.F90                        \
	basic/varinfo.F90                      \
	basic/write_iter.F90

all_c_srcs =                                   \
	basic/cublas.cc                        \
	basic/cuda_low.cc                      \
	basic/getopt_f.c                       \
	basic/io_binary.c                      \
	basic/oct_f.c                          \
	basic/signals.c                        \
	basic/sort_low.cc                      \
	basic/varia.c                          \
	basic/varinfo_low.c                    \
	basic/write_iter_low.cc

noinst_HEADERS +=                  \
	basic/accel_inc.F90        \
	basic/comm_inc.F90         \
	basic/io_binary.h          \
        basic/io_binary_f_inc.F90  \
	basic/loct_pointer_inc.F90 \
	basic/mpi_debug_inc.F90    \
	basic/mpi_lib_inc.F90      \
	basic/sort_inc.F90         \
	basic/varia.h

# Sources in basis_set
all_f_srcs +=               \
	basis_set/atomic_orbital.F90   \
	basis_set/orbitalset.F90

noinst_HEADERS +=                          \
	basis_set/atomic_orbital_inc.F90   \
	basis_set/orbitalset_inc.F90

# Sources in math
all_f_srcs +=                                 \
	math/accel_blas.F90                   \
	math/blas.F90                         \
	math/checksum_interface.F90           \
	math/fft.F90                          \
	math/fftw.F90                         \
	math/lalg_adv.F90                     \
	math/lalg_basic.F90                   \
	math/lapack.F90                       \
	math/loct_math.F90                    \
	math/math.F90                         \
	math/metis.F90                        \
	math/pblas.F90                        \
	math/pfft.F90                         \
	math/quickrnd.F90                     \
	math/scalapack.F90                    \
	math/sgfft.F90                        \
	math/solvers.F90                      \
	math/spline_filter.F90                \
	math/splines.F90                      \
	math/symm_op.F90                      \
	math/tdfunction.F90

all_c_srcs +=                                 \
	math/checksum.c                       \
	math/cufft.cc                         \
	math/dablas.c                         \
	math/fftw_low.c                       \
	math/metis_f.c                        \
	math/oct_gsl_f.c                      \
	math/spline_low.cc                    \
	math/ylm.c

noinst_HEADERS +=                             \
	math/accel_blas_inc.F90               \
	math/fft_inc.F90                      \
	math/lalg_basic_blas_inc.F90          \
	math/lalg_adv_lapack_inc.F90          \
	math/math_inc.F90                     \
	math/solvers_inc.F90                  \
	math/symm_op_inc.F90

# Sources in species
all_f_srcs +=                                    \
	species/atomic.F90                       \
	species/element.F90                      \
	species/logrid.F90                       \
	species/ps.F90                           \
	species/ps_in_grid.F90                   \
	species/ps_psf.F90                       \
	species/ps_psf_file.F90                  \
	species/ps_xml.F90                       \
	species/pseudo.F90                       \
	species/pseudo_set.F90                   \
	species/share_directory.F90              \
	species/species.F90

all_c_srcs +=                                    \
	species/element_low.cc                   \
	species/pseudo_low.cc                    \
	species/pseudo_set_low.cc                \
	species/share_directory_low.cc           \
	species/spline.cc

noinst_HEADERS +=                                \
	species/anygrid.hpp                      \
	species/base.hpp                         \
	species/element.hpp                      \
	species/detect_format.hpp                \
	species/psml.hpp                         \
	species/psp8.hpp                         \
	species/qso.hpp                          \
	species/share_directory.hpp              \
	species/set.hpp                          \
	species/spline.h                         \
	species/upf.hpp                          \
	species/upf1.hpp                         \
	species/upf2.hpp


# Sources in ions
all_f_srcs +=                                 \
	ions/atom.F90                         \
	ions/box.F90                          \
	ions/box_union.F90                    \
	ions/geometry.F90                     \
	ions/ion_dynamics.F90                 \
	ions/kpoints.F90                      \
	ions/periodic_copy.F90                \
	ions/read_coords.F90                  \
	ions/simul_box.F90                    \
	ions/spglib_f08.f90                   \
	ions/symmetries.F90

# Sources in grid
all_f_srcs +=                                 \
	grid/batch.F90                        \
	grid/batch_ops.F90                    \
	grid/boundaries.F90                   \
	grid/cube.F90                         \
	grid/cube_function.F90                \
	grid/derivatives.F90                  \
	grid/double_grid.F90                  \
	grid/fourier_space.F90                \
	grid/grid.F90                         \
	grid/index.F90                        \
	grid/io_function.F90                  \
	grid/mesh_batch.F90                   \
	grid/mesh_cube_map.F90                \
	grid/mesh_cube_parallel_map.F90       \
	grid/mesh.F90                         \
	grid/mesh_function.F90                \
	grid/mesh_init.F90                    \
	grid/mesh_interpolation.F90           \
	grid/mesh_partition.F90               \
	grid/multigrid.F90                    \
	grid/nl_operator.F90                  \
	grid/operate_f.F90                    \
	grid/partition.F90                    \
	grid/partition_transfer.F90           \
	grid/par_vec.F90                      \
	grid/restart.F90                      \
	grid/stencil_cube.F90                 \
	grid/stencil.F90                      \
	grid/stencil_star.F90                 \
	grid/stencil_stargeneral.F90          \
	grid/stencil_starplus.F90             \
	grid/stencil_variational.F90          \
	grid/subarray.F90                     \
	grid/submesh.F90                      \
	grid/symmetrizer.F90                  \
	grid/transfer_table.F90               \
	grid/vtk.F90

all_c_srcs +=                                 \
	grid/hilbert.c                        \
	grid/operate.c

noinst_HEADERS +=                             \
	grid/batch_inc.F90                    \
	grid/batch_ops_inc.F90                \
	grid/boundaries_inc.F90               \
	grid/cube_function_inc.F90            \
	grid/derivatives_inc.F90              \
	grid/double_grid_apply_inc.F90        \
	grid/fourier_space_inc.F90            \
	grid/io_function_inc.F90              \
	grid/mesh_batch_inc.F90               \
	grid/mesh_function_inc.F90            \
	grid/mesh_interpolation_inc.F90       \
	grid/multigrid_inc.F90                \
	grid/nl_operator_inc.F90              \
	grid/operate_inc.c                    \
	grid/partition_transfer_inc.F90       \
	grid/par_vec_inc.F90                  \
	grid/restart_inc.F90                  \
	grid/subarray_inc.F90                 \
	grid/submesh_inc.F90                  \
	grid/symmetrizer_inc.F90              \
	grid/vtk_inc.F90


# Sources in poisson
all_f_srcs +=                                    \
	poisson/poisson_cg.F90                   \
	poisson/poisson_corrections.F90          \
	poisson/poisson_cutoff.F90               \
	poisson/poisson.F90                      \
	poisson/poisson_fft.F90                  \
	poisson/poisson_isf.F90                  \
	poisson/poisson_multigrid.F90            \
	poisson/scaling_function.F90

all_c_srcs +=                                    \
	poisson/poisson_cutoffs.c

noinst_HEADERS +=                                \
	poisson/lazy_8_inc.F90                   \
	poisson/poisson_inc.F90

# Sources in states
all_f_srcs +=                                   \
	states/density.F90                      \
	states/smear.F90                        \
	states/states_calc.F90                  \
	states/states_dim.F90                   \
	states/states.F90                       \
	states/states_group.F90                 \
	states/states_parallel.F90              \
	states/states_restart.F90

noinst_HEADERS +=                               \
	states/density_inc.F90                  \
	states/states_calc_inc.F90              \
	states/states_inc.F90                   \
	states/states_parallel_inc.F90

# Sources in hamiltonian
all_f_srcs +=                                        \
	hamiltonian/boundary_op.F90                  \
	hamiltonian/energy.F90                       \
	hamiltonian/epot.F90                         \
	hamiltonian/functionals.F90                  \
	hamiltonian/gauge_field.F90                  \
	hamiltonian/hamiltonian_base.F90             \
	hamiltonian/hamiltonian.F90                  \
	hamiltonian/ion_interaction.F90              \
	hamiltonian/kb_projector.F90                 \
	hamiltonian/kick.F90                         \
	hamiltonian/projector.F90                    \
	hamiltonian/projector_matrix.F90             \
	hamiltonian/species_pot.F90                  \
	hamiltonian/xc.F90

noinst_HEADERS +=                                    \
        hamiltonian/functionals_list.F90             \
	hamiltonian/hamiltonian_base_inc.F90         \
	hamiltonian/hamiltonian_inc.F90              \
	hamiltonian/kb_projector_inc.F90             \
	hamiltonian/projector_inc.F90

# Sources in system
all_f_srcs +=                                   \
	system/current.F90                      \
	system/eigen_cg.F90                     \
	system/eigen_rmmdiis.F90                \
	system/eigensolver.F90                  \
	system/energy_calc.F90                  \
	system/exponential.F90                  \
	system/forces.F90                       \
	system/output.F90                       \
	system/preconditioners.F90              \
	system/subspace.F90                     \
	system/system.F90                       \
	system/v_ks.F90

noinst_HEADERS +=                               \
	system/eigen_cg_inc.F90                 \
	system/eigen_evolution_inc.F90          \
	system/eigen_plan_inc.F90               \
	system/eigen_rmmdiis_inc.F90            \
	system/energy_calc_inc.F90              \
	system/forces_inc.F90                   \
	system/output_h_inc.F90                 \
	system/output_states_inc.F90            \
	system/preconditioners_inc.F90          \
	system/subspace_inc.F90

# Sources in scf
all_f_srcs +=                                \
	scf/ground_state.F90                 \
	scf/lcao.F90                         \
	scf/mix.F90                          \
	scf/mixing_metric.F90                \
	scf/scf.F90

noinst_HEADERS +=                            \
	scf/lcao_inc.F90                     \
	scf/mix_inc.F90


# Sources in td
all_f_srcs +=                               \
	td/potential_interpolation.F90      \
	td/propagator_base.F90              \
	td/propagator_etrs.F90              \
	td/propagator_expmid.F90            \
	td/propagator.F90                   \
	td/td.F90                           \
	td/td_write.F90

all_srcs = $(all_f_srcs) $(all_c_srcs)

include/options.h: $(all_f_srcs) $(noinst_HEADERS)
	$(top_srcdir)/build/mk_varinfo.pl -s $(top_srcdir) -b $(top_builddir)

include/defaults.h: $(all_f_srcs) $(noinst_HEADERS)
	$(top_srcdir)/build/mk_varinfo.pl -s $(top_srcdir) -b $(top_builddir)

bin_PROGRAMS = pulpito

pulpito_f_srcs =                             \
	main/test.F90                        \
	main/run.F90                         \
	main/main.F90

pulpito_SOURCES = $(pulpito_f_srcs) $(all_srcs)
pulpito_LDADD = $(all_LIBS)
pulpito_LINK = $(FCLINK)


CLEANFILES += pulpito

EXEEXT =

EXTRA_PROGRAMS = pulpito-help

bin_PROGRAMS += @UTILITY_PROGRAMS@

pulpito_help_SOURCES = $(all_srcs) main/help.F90
pulpito_help_LINK = $(FCLINK)
pulpito_help_LDADD = $(all_LIBS)

# debug files from MacOS gcc with -g
clean-local: clean-fdep
	rm -rf *.dSYM
	rm -f include/config_F90.h include/options.h include/defaults.h


FORTRAN_CPP=@FCCPP@ @CPPFLAGS@ $(AM_CPPFLAGS)
PERL_ARGS="@DEBUG@" "@F90_ACCEPTS_LINE_NUMBERS@" "@F90_FORALL@"
HEADERS_DEP=include/config_F90.h include/options.h include/defaults.h

EXTRA_DIST = fdep/fortran_dependencies.mk \
  fdep/fortran_dependencies.pl

@FORTRAN_MODULE_DEPS@
