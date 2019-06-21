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
basic_f_srcs =                                 \
	basic/accel.F90                        \
	basic/alloc_cache.F90                  \
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
	basic/iihash.F90                       \
	basic/io_binary_f.F90                  \
	basic/io_csv_f.F90                     \
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
	basic/walltimer.F90                    \
	basic/write_iter.F90

basic_c_srcs =                                 \
	basic/alloc_cache_low.cc               \
	basic/cublas.cc                        \
	basic/cuda_low.cc                      \
	basic/getopt_f.c                       \
	basic/iihash_low.cc                    \
	basic/io_binary.c                      \
	basic/io_csv.c                         \
	basic/oct_f.c                          \
	basic/oct_gdlib_f.c                    \
	basic/recipes.c                        \
	basic/signals.c                        \
	basic/sort_low.cc                      \
	basic/varia.c                          \
	basic/varinfo_low.c                    \
	basic/write_iter_low.cc

basic_srcs = $(basic_f_srcs) $(basic_c_srcs)

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
basis_set_f_srcs =                     \
	basis_set/atomic_orbital.F90   \
	basis_set/loewdin.F90          \
	basis_set/orbitalbasis.F90     \
	basis_set/orbitalset.F90       \
	basis_set/orbitalset_utils.F90

basis_set_srcs = $(basis_set_f_srcs)

noinst_HEADERS +=                          \
	basis_set/atomic_orbital_inc.F90   \
	basis_set/loewdin_inc.F90          \
	basis_set/orbitalbasis_inc.F90     \
	basis_set/orbitalset_inc.F90       \
	basis_set/orbitalset_utils_inc.F90


# Sources in math
math_f_srcs =                                 \
	math/accel_blas.F90                   \
	math/blas.F90                         \
	math/checksum_interface.F90           \
	math/compressed_sensing.F90           \
	math/fft.F90                          \
	math/fftw.F90                         \
	math/gauss_legendre.F90               \
	math/lalg_adv.F90                     \
	math/lalg_basic.F90                   \
	math/lapack.F90                       \
	math/loct_math.F90                    \
	math/math.F90                         \
	math/matrix.F90                       \
	math/matrix_ops.F90                   \
	math/metis.F90                        \
	math/minimizer.F90                    \
	math/nfft.F90                         \
	math/ode_solver.F90                   \
	math/pblas.F90                        \
	math/permutations.F90                 \
	math/pfft.F90                         \
	math/pnfft.F90                        \
	math/qshep.F90                        \
	math/quickrnd.F90                     \
	math/root_solver.F90                  \
	math/scalapack.F90                    \
	math/sgfft.F90                        \
	math/solvers.F90                      \
	math/sparskit.F90                     \
	math/spline_filter.F90                \
	math/splines.F90                      \
	math/symm_op.F90                      \
	math/tdfunction.F90                   \
	math/young.F90

math_c_srcs =                                 \
	math/checksum.c                       \
	math/cufft.cc                         \
	math/dablas.c                         \
	math/fftw_low.c                       \
	math/metis_f.c                        \
	math/minimizer_low.c                  \
	math/nfft_f.c                         \
	math/oct_gsl_f.c                      \
	math/spline_low.cc                    \
	math/ylm.c

math_srcs = $(math_f_srcs) $(math_c_srcs)

noinst_HEADERS +=                             \
	math/accel_blas_inc.F90               \
	math/fft_inc.F90                      \
	math/lalg_basic_blas_inc.F90          \
	math/lalg_adv_lapack_inc.F90          \
	math/math_inc.F90                     \
	math/matrix_inc.F90                   \
	math/matrix_ops_inc.F90               \
	math/nfft_inc.F90                     \
	math/ode_solver_inc.F90               \
	math/pnfft_inc.F90                    \
	math/root_solver_inc.F90              \
	math/sparskit_inc.F90                 \
	math/solvers_inc.F90                  \
	math/symm_op_inc.F90

# Sources in species
species_f_srcs =                                 \
	species/atomic.F90                       \
	species/element.F90                      \
	species/logrid.F90                       \
	species/ps_cpi.F90                       \
	species/ps_cpi_file.F90                  \
	species/ps.F90                           \
	species/ps_fhi.F90                       \
	species/ps_fhi_file.F90                  \
	species/ps_hgh.F90                       \
	species/ps_in_grid.F90                   \
	species/ps_psf.F90                       \
	species/ps_psf_file.F90                  \
	species/ps_xml.F90                       \
	species/pseudo.F90                       \
	species/pseudo_set.F90                   \
	species/share_directory.F90              \
	species/species.F90

species_c_srcs =                                 \
	species/element_low.cc                   \
	species/pseudo_low.cc                    \
	species/pseudo_set_low.cc                \
	species/share_directory_low.cc           \
	species/spline.cc

species_srcs = $(species_f_srcs) $(species_c_srcs)

noinst_HEADERS +=                                \
	species/anygrid.hpp                      \
	species/base.hpp                         \
	species/element.hpp                      \
	species/detect_format.hpp                \
	species/ps_pspio_inc.F90                 \
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
ions_f_srcs =                                 \
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
	ions/symmetries.F90                   \
	ions/vibrations.F90                   \
	ions/xyz_adjust.F90

ions_c_srcs =                                 \
	ions/symmetries_finite.c

ions_srcs = $(ions_f_srcs) $(ions_c_srcs)


# Sources in grid
grid_f_srcs =                                 \
	grid/basins.F90                       \
	grid/batch.F90                        \
	grid/batch_ops.F90                    \
	grid/boundaries.F90                   \
	grid/cube.F90                         \
	grid/cube_function.F90                \
	grid/curv_briggs.F90                  \
	grid/curv_gygi.F90                    \
	grid/curvilinear.F90                  \
	grid/curv_modine.F90                  \
	grid/derivatives.F90                  \
	grid/double_grid.F90                  \
	grid/fourier_shell.F90                \
	grid/fourier_space.F90                \
	grid/grid.F90                         \
	grid/hypercube.F90                    \
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
	grid/modelmb_1part.F90                \
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
	grid/volume.F90                       \
	grid/vtk.F90

grid_c_srcs =                                 \
	grid/hilbert.c                        \
	grid/operate.c

grid_srcs = $(grid_f_srcs) $(grid_c_srcs)

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
poisson_f_srcs =                                 \
	poisson/poisson_cg.F90                   \
	poisson/poisson_corrections.F90          \
	poisson/poisson_cutoff.F90               \
	poisson/poisson.F90                      \
	poisson/poisson_fft.F90                  \
	poisson/poisson_fmm.F90                  \
	poisson/poisson_isf.F90                  \
	poisson/poisson_libisf.F90               \
	poisson/poisson_multigrid.F90            \
	poisson/poisson_no.F90                   \
	poisson/scaling_function.F90

poisson_c_srcs =                                 \
	poisson/poisson_cutoffs.c

poisson_srcs = $(poisson_f_srcs) $(poisson_c_srcs)

noinst_HEADERS +=                                \
	poisson/lazy_8_inc.F90                   \
	poisson/poisson_inc.F90                  \
	poisson/poisson_direct_sm_inc.F90        \
	poisson/poisson_init_direct_inc.F90      \
	poisson/solver_1d_solve_inc.F90

# Sources in states
states_f_srcs =                                 \
	states/density.F90                      \
	states/excited_states.F90               \
	states/linear_response.F90              \
	states/modelmb_density_matrix.F90       \
	states/modelmb_exchange_syms.F90        \
	states/modelmb_particles.F90            \
	states/opt_control_state.F90            \
	states/physics_op.F90                   \
	states/smear.F90                        \
	states/states_block.F90                 \
	states/states_calc.F90                  \
	states/states_dim.F90                   \
	states/states.F90                       \
	states/states_group.F90                 \
	states/states_io.F90                    \
	states/states_parallel.F90              \
	states/states_restart.F90

states_srcs = $(states_f_srcs)

noinst_HEADERS +=                               \
	states/density_inc.F90                  \
	states/excited_states_inc.F90           \
	states/linear_response_inc.F90          \
	states/modelmb_density_matrix_inc.F90   \
	states/modelmb_exchange_syms_inc.F90    \
	states/physics_op_inc.F90               \
	states/states_block_inc.F90             \
	states/states_calc_inc.F90              \
	states/states_inc.F90                   \
	states/states_parallel_inc.F90

# Sources in hamiltonian
hamiltonian_f_srcs =                                 \
	hamiltonian/berry.F90                        \
	hamiltonian/born_charges.F90                 \
	hamiltonian/boundary_op.F90                  \
	hamiltonian/energy.F90                       \
	hamiltonian/epot.F90                         \
	hamiltonian/functionals.F90                  \
	hamiltonian/gauge_field.F90                  \
	hamiltonian/hamiltonian_base.F90             \
	hamiltonian/hamiltonian.F90                  \
	hamiltonian/hgh_projector.F90                \
	hamiltonian/hirshfeld.F90                    \
	hamiltonian/ion_interaction.F90              \
	hamiltonian/kb_projector.F90                 \
	hamiltonian/kick.F90                         \
	hamiltonian/lasers.F90                       \
	hamiltonian/lda_u.F90                        \
	hamiltonian/lda_u_io.F90                     \
	hamiltonian/libvdwxc.F90                     \
	hamiltonian/oct_exchange.F90                 \
	hamiltonian/pcm.F90                          \
	hamiltonian/pcm_eom.F90                      \
	hamiltonian/projector.F90                    \
	hamiltonian/projector_matrix.F90             \
	hamiltonian/rkb_projector.F90                \
	hamiltonian/scdm.F90                         \
	hamiltonian/scf_tol.F90                      \
	hamiltonian/scissor.F90                      \
	hamiltonian/species_pot.F90                  \
	hamiltonian/vdw_ts.F90                       \
	hamiltonian/xc.F90

hamiltonian_c_srcs =                                 \
	hamiltonian/vdw_ts_low.c

hamiltonian_srcs = $(hamiltonian_f_srcs) $(hamiltonian_c_srcs)

noinst_HEADERS +=                                    \
        hamiltonian/dft_u_noncollinear_inc.F90       \
        hamiltonian/functionals_list.F90             \
        hamiltonian/fxc_inc.F90                      \
	hamiltonian/hamiltonian_base_inc.F90         \
	hamiltonian/hamiltonian_inc.F90              \
	hamiltonian/hgh_projector_inc.F90            \
	hamiltonian/kb_projector_inc.F90             \
	hamiltonian/kxc_inc.F90                      \
	hamiltonian/lasers_inc.F90                   \
	hamiltonian/lda_u_inc.F90                    \
	hamiltonian/oct_exchange_inc.F90             \
	hamiltonian/projector_inc.F90                \
	hamiltonian/scdm_inc.F90                     \
	hamiltonian/scissor_inc.F90                  \
	hamiltonian/vxc_inc.F90

# Sources in system
system_f_srcs =                                 \
	system/current.F90                      \
	system/dos.F90                          \
	system/eigen_cg.F90                     \
	system/eigen_lobpcg.F90                 \
	system/eigen_rmmdiis.F90                \
	system/eigensolver.F90                  \
	system/elf.F90                          \
	system/energy_calc.F90                  \
	system/exponential.F90                  \
	system/forces.F90                       \
	system/linear_solver.F90                \
	system/magnetic.F90                     \
	system/output.F90                       \
	system/output_me.F90                    \
	system/partial_charges.F90              \
	system/pert.F90                         \
	system/preconditioners.F90              \
	system/stress.F90                       \
	system/subspace.F90                     \
	system/system.F90                       \
	system/v_ks.F90                         \
	system/xc_ks_inversion.F90              \
	system/xc_oep.F90

system_srcs = $(system_f_srcs)

noinst_HEADERS +=                               \
	system/eigen_cg_inc.F90                 \
	system/eigen_evolution_inc.F90          \
	system/eigen_lobpcg_inc.F90             \
	system/eigen_plan_inc.F90               \
	system/eigen_rmmdiis_inc.F90            \
	system/energy_calc_inc.F90              \
	system/forces_inc.F90                   \
	system/linear_solver_inc.F90            \
	system/output_berkeleygw_inc.F90        \
	system/output_etsf_inc.F90              \
	system/output_h_inc.F90                 \
	system/output_linear_response_inc.F90   \
	system/output_me_inc.F90                \
	system/output_modelmb_inc.F90           \
	system/output_states_inc.F90            \
	system/pert_inc.F90                     \
	system/preconditioners_inc.F90          \
	system/subspace_inc.F90                 \
	system/xc_kli_inc.F90                   \
	system/xc_kli_pauli_inc.F90             \
	system/xc_oep_inc.F90                   \
	system/xc_oep_sic_inc.F90               \
	system/xc_oep_x_inc.F90

# Sources in scf
scf_f_srcs =                                 \
	scf/ground_state.F90                 \
	scf/lcao.F90                         \
	scf/lda_u_mixer.F90                  \
	scf/mix.F90                          \
	scf/mixing_metric.F90                \
	scf/rdmft.F90                        \
	scf/scf.F90

scf_srcs = $(scf_f_srcs)

noinst_HEADERS +=                            \
	scf/lcao_inc.F90                     \
	scf/mix_inc.F90


# Sources in td
td_f_srcs =                                 \
	td/pes.F90                          \
	td/pes_flux.F90                     \
	td/pes_mask.F90                     \
	td/pes_out.F90                      \
	td/pes_spm.F90                      \
	td/potential_interpolation.F90      \
	td/propagator_base.F90              \
	td/propagator_cn.F90                \
	td/propagator_etrs.F90              \
	td/propagator_expmid.F90            \
	td/propagator.F90                   \
	td/propagator_magnus.F90            \
	td/propagator_qoct.F90              \
	td/propagator_rk.F90                \
	td/spectrum.F90                     \
	td/td_calc.F90                      \
	td/td.F90                           \
	td/td_write.F90

td_srcs = $(td_f_srcs)

noinst_HEADERS +=                           \
	td/pes_mask_out_inc.F90             \
	td/pes_flux_out_inc.F90


# Sources in opt_control
opt_control_f_srcs =                                 \
	opt_control/controlfunction.F90              \
	opt_control/filter.F90                       \
	opt_control/initst.F90                       \
	opt_control/opt_control.F90                  \
	opt_control/opt_control_global.F90           \
	opt_control/opt_control_iter.F90             \
	opt_control/propagation.F90                  \
	opt_control/target.F90

opt_control_srcs = $(opt_control_f_srcs)

noinst_HEADERS +=                                    \
	opt_control/check_input_inc.F90              \
	opt_control/controlfunction_trans_inc.F90    \
	opt_control/finalcheck_inc.F90               \
	opt_control/opt_control_c_inc.F90            \
	opt_control/target_classical_inc.F90         \
	opt_control/target_density_inc.F90           \
	opt_control/target_excited_inc.F90           \
	opt_control/target_exclude_inc.F90           \
	opt_control/target_groundstate_inc.F90       \
	opt_control/target_gstransformation_inc.F90  \
	opt_control/target_hhg_inc.F90               \
	opt_control/target_local_inc.F90             \
	opt_control/target_spin_inc.F90              \
	opt_control/target_tdlocal_inc.F90           \
	opt_control/target_userdefined_inc.F90       \
	opt_control/target_velocity_inc.F90


# Sources in sternheimer
sternheimer_f_srcs =                                 \
	sternheimer/em_resp_calc.F90                 \
	sternheimer/em_resp.F90                      \
	sternheimer/kdotp_calc.F90                   \
	sternheimer/kdotp.F90                        \
	sternheimer/phonons_lr.F90                   \
	sternheimer/sternheimer.F90                  \
	sternheimer/vdw.F90

sternheimer_srcs = $(sternheimer_f_srcs)

noinst_HEADERS +=                                    \
	sternheimer/em_resp_calc_inc.F90             \
	sternheimer/em_resp_inc.F90                  \
	sternheimer/kdotp_calc_inc.F90               \
	sternheimer/phonons_lr_inc.F90               \
	sternheimer/sternheimer_inc.F90


all_f_srcs = $(basic_f_srcs) $(basis_set_f_srcs) $(math_f_srcs) $(species_f_srcs) \
	     $(ions_f_srcs) $(grid_f_srcs) $(poisson_f_srcs) $(states_f_srcs) \
	     $(hamiltonian_f_srcs) $(system_f_srcs) $(scf_f_srcs) $(td_f_srcs) \
	     $(opt_control_f_srcs) $(sternheimer_f_srcs)

all_c_srcs = $(basic_c_srcs) $(math_c_srcs) $(species_c_srcs) $(ions_c_srcs) \
	     $(grid_c_srcs) $(poisson_c_srcs) $(hamiltonian_c_srcs)

all_srcs = $(all_f_srcs) $(all_c_srcs)

include/options.h: $(all_f_srcs) $(noinst_HEADERS)
	$(top_srcdir)/build/mk_varinfo.pl -s $(top_srcdir) -b $(top_builddir)

include/defaults.h: $(all_f_srcs) $(noinst_HEADERS)
	$(top_srcdir)/build/mk_varinfo.pl -s $(top_srcdir) -b $(top_builddir)

bin_PROGRAMS = octopus

octopus_f_srcs =                             \
	main/invert_ks.F90                   \
	main/memory.F90                      \
	main/unocc.F90                       \
	main/static_pol.F90                  \
	main/geom_opt.F90                    \
	main/phonons_fd.F90                  \
	main/casida.F90                      \
	main/pulpo.F90                       \
	main/test.F90                        \
	main/run.F90                         \
	main/main.F90

noinst_HEADERS += \
	main/casida_inc.F90

octopus_SOURCES = $(octopus_f_srcs) $(all_srcs)
octopus_LDADD = $(all_LIBS)
octopus_LINK = $(FCLINK)


CLEANFILES += octopus

EXEEXT =

EXTRA_PROGRAMS = oct-propagation_spectrum \
	         oct-harmonic-spectrum \
                 oct-casida_spectrum      \
                 oct-vibrational_spectrum \
                 oct-infrared_spectrum    \
                 oct-xyz-anim \
                 oct-center-geom \
                 oct-help \
                 oct-oscillator-strength \
                 oct-dielectric-function \
                 oct-photoelectron_spectrum \
                 oct-convert \
                 oct-local_multipoles \
                 oct-conductivity \
                 oct-floquet \
                 oct-unfold

bin_PROGRAMS += @UTILITY_PROGRAMS@

oct_propagation_spectrum_SOURCES = $(td_srcs) $(scf_srcs) $(system_srcs) $(hamiltonian_srcs) \
				$(basis_set_srcs) $(poisson_srcs) $(states_srcs) $(grid_srcs) \
				$(ions_srcs) $(species_srcs) $(math_srcs) $(basic_srcs) \
				utils/propagation_spectrum.F90
oct_propagation_spectrum_LINK = $(FCLINK)
oct_propagation_spectrum_LDADD = $(all_LIBS)
oct_oscillator_strength_SOURCES = $(td_srcs) $(scf_srcs) $(system_srcs) $(hamiltonian_srcs) \
				$(basis_set_srcs) $(poisson_srcs) $(states_srcs) $(grid_srcs) \
				$(ions_srcs) $(species_srcs) $(math_srcs) $(basic_srcs) \
				utils/oscillator_strength.F90
oct_oscillator_strength_LINK = $(FCLINK)
oct_oscillator_strength_LDADD = $(all_LIBS)
oct_dielectric_function_SOURCES = $(td_srcs) $(scf_srcs) $(system_srcs) $(hamiltonian_srcs) \
				$(basis_set_srcs) $(poisson_srcs) $(states_srcs) $(grid_srcs) \
				$(ions_srcs) $(species_srcs) $(math_srcs) $(basic_srcs) \
				utils/dielectric_function.F90
oct_dielectric_function_LINK = $(FCLINK)
oct_dielectric_function_LDADD = $(all_LIBS)
oct_harmonic_spectrum_SOURCES = $(td_srcs) $(scf_srcs) $(system_srcs) $(hamiltonian_srcs) \
				$(basis_set_srcs) $(poisson_srcs) $(states_srcs) $(grid_srcs) \
				$(ions_srcs) $(species_srcs) $(math_srcs) $(basic_srcs) \
				utils/harmonic_spectrum.F90
oct_harmonic_spectrum_LINK = $(FCLINK)
oct_harmonic_spectrum_LDADD = $(all_LIBS)
oct_photoelectron_spectrum_SOURCES = $(td_srcs) $(scf_srcs) $(system_srcs) $(hamiltonian_srcs) \
				$(basis_set_srcs) $(poisson_srcs) $(states_srcs) $(grid_srcs) \
				$(ions_srcs) $(species_srcs) $(math_srcs) $(basic_srcs) \
				utils/photoelectron_spectrum.F90
oct_photoelectron_spectrum_LINK = $(FCLINK)
oct_photoelectron_spectrum_LDADD = $(all_LIBS)
oct_xyz_anim_SOURCES = $(grid_srcs) $(ions_srcs) $(species_srcs) $(math_srcs) $(basic_srcs) \
		utils/xyzanim.F90
oct_xyz_anim_LINK = $(FCLINK)
oct_xyz_anim_LDADD = $(all_LIBS)
oct_casida_spectrum_SOURCES = $(grid_srcs) $(ions_srcs) $(species_srcs) $(math_srcs) $(basic_srcs) \
			utils/casida_spectrum.F90
oct_casida_spectrum_LINK = $(FCLINK)
oct_casida_spectrum_LDADD = $(all_LIBS)
oct_vibrational_spectrum_SOURCES = $(td_srcs) $(scf_srcs) $(system_srcs) $(hamiltonian_srcs) \
			$(basis_set_srcs) $(poisson_srcs) $(states_srcs) $(grid_srcs) \
			$(ions_srcs) $(species_srcs) $(math_srcs) $(basic_srcs) \
			utils/vibrational.F90
oct_vibrational_spectrum_LINK = $(FCLINK)
oct_vibrational_spectrum_LDADD = $(all_LIBS)
oct_infrared_spectrum_SOURCES = $(grid_srcs) $(ions_srcs) $(species_srcs) $(math_srcs) $(basic_srcs) \
				utils/infrared.F90
oct_infrared_spectrum_LINK = $(FCLINK)
oct_infrared_spectrum_LDADD = $(all_LIBS)
oct_center_geom_SOURCES = $(grid_srcs) $(ions_srcs) $(species_srcs) $(math_srcs) $(basic_srcs) \
			  utils/centergeom.F90
oct_center_geom_LINK = $(FCLINK)
oct_center_geom_LDADD = $(all_LIBS)
oct_convert_SOURCES = $(td_srcs) $(scf_srcs) $(system_srcs) $(hamiltonian_srcs) \
			$(basis_set_srcs) $(poisson_srcs) $(states_srcs) $(grid_srcs) \
			$(ions_srcs) $(species_srcs) $(math_srcs) $(basic_srcs) \
			utils/convert.F90
oct_convert_LINK = $(FCLINK)
oct_convert_LDADD = $(all_LIBS)
oct_local_multipoles_SOURCES = $(td_srcs) $(scf_srcs) $(system_srcs) $(hamiltonian_srcs) \
			$(basis_set_srcs) $(poisson_srcs) $(states_srcs) $(grid_srcs) \
			$(ions_srcs) $(species_srcs) $(math_srcs) $(basic_srcs) \
			utils/local_multipoles.F90 \
			utils/local_write.F90
oct_local_multipoles_LINK = $(FCLINK)
oct_local_multipoles_LDADD = $(all_LIBS)
oct_conductivity_SOURCES = $(td_srcs) $(scf_srcs) $(system_srcs) $(hamiltonian_srcs) \
			$(basis_set_srcs) $(poisson_srcs) $(states_srcs) $(grid_srcs) \
			$(ions_srcs) $(species_srcs) $(math_srcs) $(basic_srcs) \
			utils/conductivity.F90
oct_conductivity_LINK = $(FCLINK)
oct_conductivity_LDADD = $(all_LIBS)
oct_help_SOURCES = $(basic_srcs) utils/help.F90
oct_help_LINK = $(FCLINK)
oct_help_LDADD = $(all_LIBS)
oct_floquet_SOURCES = $(td_srcs) $(scf_srcs) $(system_srcs) $(hamiltonian_srcs) $(basis_set_srcs) \
			$(poisson_srcs) $(states_srcs) $(grid_srcs) \
			$(ions_srcs) $(species_srcs) $(math_srcs) $(basic_srcs) \
			utils/floquet.F90
oct_floquet_LINK = $(FCLINK)
oct_floquet_LDADD = $(all_LIBS)

oct_unfold_SOURCES    = $(td_srcs) $(scf_srcs) $(system_srcs) $(hamiltonian_srcs) \
                        $(basis_set_srcs) $(poisson_srcs) $(states_srcs) $(grid_srcs) \
                        $(ions_srcs) $(species_srcs) $(math_srcs) $(basic_srcs) \
                        utils/unfold.F90
oct_unfold_LINK = $(FCLINK)
oct_unfold_LDADD = $(all_LIBS)



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
