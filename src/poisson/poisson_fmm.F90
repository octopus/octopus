!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch, 
!! J. Alberdi-Rodriguez, P. Garcia Risueño, M. Oliveira
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"
#ifdef HAVE_LIBFM
#include "fcs_fconfig.h"
#endif

module poisson_fmm_oct_m
  use boundaries_oct_m
  use box_parallelepiped_oct_m
  use cube_oct_m
  use derivatives_oct_m
  use fft_oct_m
  use global_oct_m
  use index_oct_m
  use io_function_oct_m
  use iso_c_binding
  use lalg_basic_oct_m
  use loct_math_oct_m
  use math_oct_m
  use mesh_oct_m
  use mesh_cube_parallel_map_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use multicomm_oct_m
  use nl_operator_oct_m
  use par_vec_oct_m
  use parser_oct_m
  use profiling_oct_m
  use simul_box_oct_m
  use space_oct_m
  use stencil_star_oct_m

#ifdef HAVE_LIBFM
  use fcs_module
  use iso_fortran_env
#endif

  implicit none

  private

  public ::                      &
    poisson_fmm_t,               &
    poisson_fmm_init,            &
    poisson_fmm_end,             &
    poisson_fmm_solve

#ifndef HAVE_LIBFM
  integer, parameter :: fcs_integer_kind_isoc = 4
#endif

  type poisson_fmm_t
    private
    FLOAT   :: delta_E_fmm
    FLOAT   :: alpha_fmm  !< Alpha for the correction of the FMM
    type(mpi_grp_t) :: all_nodes_grp !< The communicator for all nodes.
    type(mpi_grp_t) :: perp_grp      !< The communicator perpendicular to the mesh communicator.
    integer(kind = fcs_integer_kind_isoc) :: nlocalcharges
    integer    :: sp !< Local start point
    integer    :: ep !< Local end point
    integer, allocatable :: disps(:)
    integer, allocatable :: dsize(:) !< Local size
    type(nl_operator_t) :: corrector
    type(derivatives_t), pointer :: der
    type(c_ptr) ::  handle !< The FMM identifier
  end type poisson_fmm_t

contains

  ! ---------------------------------------------------------
  !> Initialises the FMM parameters and vectors. Also it calls to 
  !! the library initialisation.
  subroutine poisson_fmm_init(this, space, der, all_nodes_comm)
    type(poisson_fmm_t),         intent(out)   :: this
    type(space_t),               intent(in)    :: space
    type(derivatives_t), target, intent(in)    :: der
    integer,                     intent(in)    :: all_nodes_comm

#ifdef HAVE_LIBFM    
    integer :: is
    logical, allocatable :: remains(:)
    integer, allocatable :: dend(:)
    integer :: subcomm, cdim
    character(len = 8) ::  method
    type(c_ptr) ::  ret
    type(mesh_t), pointer :: mesh

    logical ::  short_range_flag = .true.
    real(kind = fcs_real_kind_isoc), dimension(3)   ::  box_a 
    real(kind = fcs_real_kind_isoc), dimension(3)   ::  box_b 
    real(kind = fcs_real_kind_isoc), dimension(3)   ::  box_c 
    real(kind = fcs_real_kind_isoc), dimension(3)   ::  offset = (/M_ZERO,M_ZERO,M_ZERO/)
    logical, dimension(3)                           ::  periodicity = (/.false.,.false.,.false./)
    integer(kind = fcs_integer_kind_isoc)           ::  total_particles

    integer(kind = fcs_integer_kind_isoc)           ::  local_particle_count = -1
    real(kind = fcs_real_kind_isoc), dimension(8)   ::  local_charges
    real(kind = fcs_real_kind_isoc), dimension(24)  ::  local_coordinates

    PUSH_SUB(poisson_fmm_init)

    method = "fmm"
    !%Variable DeltaEFMM
    !%Type float
    !%Default 0.0001 
    !%Section Hamiltonian::Poisson 
    !%Description
    !% Dimensionless parameter for relative convergence of <tt>PoissonSolver = FMM</tt>.
    !% Sets energy error bound.
    !% Strong inhomogeneous systems may violate the error bound.
    !% For inhomogeneous systems we have an error-controlled sequential version available
    !% (from Ivo Kabadshow).
    !%
    !% Our implementation of FMM (based on H. Dachsel, <i>J. Chem. Phys.</i> <b>131</b>,
    !% 244102 (2009)) can keep the error of the Hartree energy below an
    !% arbitrary bound. The quotient of the value chosen for the maximum
    !% error in the Hartree energy and the value of the Hartree energy is
    !% <tt>DeltaEFMM</tt>.
    !%
    !%End
    call parse_variable(parser, 'DeltaEFMM', CNST(1e-4), this%delta_E_fmm)

    !%Variable AlphaFMM
    !%Type float
    !%Default 0.291262136
    !%Section Hamiltonian::Poisson 
    !%Description
    !% Dimensionless parameter for the correction of the self-interaction of the
    !% electrostatic Hartree potential, when using <tt>PoissonSolver = FMM</tt>.
    !% 
    !% Octopus represents charge density on a real-space grid, each
    !% point containing a value <math>\rho</math> corresponding to the charge
    !% density in the cell centered in such point. Therefore, the
    !% integral for the Hartree potential at point <math>i</math>, <math>V_H(i)</math>, can be reduced to a summation:
    !%
    !% <math>V_H(i) = \frac{\Omega}{4\pi\varepsilon_0} \sum_{i \neq j}
    !% \frac{\rho(\vec{r}(j))}{|\vec{r}(j) - \vec{r}(i)|} + V_{self.int.}(i)</math>
    !% where <math>\Omega</math> is the volume element of the mesh, and <math>\vec{r}(j)</math> is the
    !% position of the point <math>j</math>. The <math>V_{self.int.}(i)</math> corresponds to
    !% the integral over the cell centered on the point <math>i</math> that is necessary to
    !% calculate the Hartree potential at point <math>i</math>:
    !%
    !% <math>V_{self.int.}(i)=\frac{1}{4\pi\varepsilon_0}
    !% \int_{\Omega(i)}d\vec{r} \frac{\rho(\vec{r}(i))}{|\vec{r}-\vec{r}(i)|}</math>
    !%
    !% In the FMM version implemented into Octopus, a correction method
    !% for <math>V_H(i)</math> is used
    !% (see Garc&iacute;a-Risue&ntilde;o <i>et al.</i>, <i>J. Comp. Chem.</i> <b>35</b>, 427 (2014)).
    !% This method defines cells neighbouring cell <math>i</math>, which
    !% have volume <math>\Omega(i)/8</math> (in 3D) and charge density obtained by
    !% interpolation. In the calculation of <math>V_H(i)</math>, in order to avoid
    !% double counting of charge, and to cancel part of the errors arising
    !% from considering the distances constant in the summation above, a
    !% term <math>-\alpha_{FMM}V_{self.int.}(i)</math> is added to the summation (see
    !% the paper for the explicit formulae).
    !%End
    call parse_variable(parser, 'AlphaFMM', CNST(0.291262136), this%alpha_fmm)

    ! FMM: Variable periodic sets periodicity
    ! 0 = open system
    ! 1 = 1D periodic system
    ! 2 = 2D periodic system
    ! 3 = 3D periodic system

    call mpi_grp_init(this%all_nodes_grp, all_nodes_comm)

    this%der => der

    if (mpi_world%size == 1) then
      cdim = 1

      SAFE_ALLOCATE(this%disps(1:1))
      SAFE_ALLOCATE(dend(1:1))
      SAFE_ALLOCATE(this%dsize(1:1))

      dend = der%mesh%np
      this%sp = 1 
      this%ep = der%mesh%np
      this%dsize(1) = der%mesh%np
      this%disps = 0
      this%nlocalcharges = this%dsize(1)

      subcomm = this%all_nodes_grp%comm

    else 
#ifdef HAVE_MPI
      call MPI_Cartdim_get(this%all_nodes_grp%comm, cdim, mpi_err)

      SAFE_ALLOCATE(remains(1:cdim))

      remains = .true.
      remains(1) = .false.

      call MPI_Cart_sub(this%all_nodes_grp%comm, remains(1), subcomm, mpi_err)
      call mpi_grp_init(this%perp_grp, subcomm)

      SAFE_ALLOCATE(this%disps(1:this%perp_grp%size))
      SAFE_ALLOCATE(dend(1:this%perp_grp%size))
      SAFE_ALLOCATE(this%dsize(1:this%perp_grp%size))

      call multicomm_divide_range(der%mesh%np, this%perp_grp%size, this%disps, dend, this%dsize)

      this%sp = this%disps(this%perp_grp%rank + 1)
      this%ep = dend(this%perp_grp%rank + 1)
      this%nlocalcharges = this%dsize(this%perp_grp%rank + 1)
      this%nlocalcharges = der%mesh%np
      this%disps = this%disps - 1

#endif
    end if

    mesh => der%mesh

    if (space%is_periodic()) then
      select type (box => mesh%sb%box)
      type is (box_parallelepiped_t)
        if (.not. all(box%half_length == box%half_length(1))) then
          call messages_not_implemented("FMM Poisson solver with non-cubic boxes.")
        end if
      class default
        call messages_not_implemented("FMM Poisson solver with non-cubic boxes.")
      end select
    end if

    total_particles = mesh%np_global
    ret = fcs_init(this%handle, trim(adjustl(method)) // c_null_char, subcomm)

    box_a = (/mesh%sb%lsize(1)*2+1,M_ZERO,M_ZERO/)
    box_b = (/M_ZERO,mesh%sb%lsize(2)*2+1,M_ZERO/)
    box_c = (/M_ZERO,M_ZERO,mesh%sb%lsize(3)*2+1/)

    ret = fcs_common_set(this%handle, short_range_flag, box_a, box_b, box_c, offset, periodicity, total_particles)

    ! this is how you set a relative error in scafacos for the FMM
    ret = fcs_fmm_set_tolerance_energy( this%handle, this%delta_E_fmm );

    call nl_operator_init(this%corrector, "FMM Correction")
    call stencil_star_get_lapl(this%corrector%stencil, der%dim, 2)
    call nl_operator_build(this%der%mesh, this%corrector, der%mesh%np, const_w = .not. this%der%mesh%use_curvilinear)

    do is = 1, this%corrector%stencil%size

      select case(sum(abs(this%corrector%stencil%points(1:der%dim, is))))
      case(0)
        this%corrector%w_re(is, 1) = CNST(27.0)/CNST(32.0) + &
          (M_ONE - this%alpha_fmm)*M_TWO*M_PI*(CNST(3.0)/(M_PI*CNST(4.0)))**(CNST(2.0)/CNST(3.0))
      case(1)
        this%corrector%w_re(is, 1) = CNST(0.0625)
      case(2)
        this%corrector%w_re(is, 1) = -CNST(0.0625)*CNST(0.25)
      end select

      this%corrector%w_re(is, 1) = this%corrector%w_re(is, 1)*der%mesh%spacing(1)*der%mesh%spacing(2)
    end do

    call nl_operator_output_weights(this%corrector)
    POP_SUB(poisson_fmm_init)
#endif
  end subroutine poisson_fmm_init

  ! ---------------------------------------------------------
  !> Release memory and call to end the library
  subroutine poisson_fmm_end(this)
    type(poisson_fmm_t), intent(inout) :: this

#ifdef HAVE_LIBFM
    type(c_ptr) ::  ret
    
    PUSH_SUB(poisson_fmm_end)

    call nl_operator_end(this%corrector)

    if (mpi_world%size > 1) call MPI_Comm_free(this%perp_grp%comm, mpi_err)
    SAFE_DEALLOCATE_A(this%disps)
    SAFE_DEALLOCATE_A(this%dsize)
    ret = fcs_destroy(this%handle)

    POP_SUB(poisson_fmm_end)
#endif
  end subroutine poisson_fmm_end

  ! ---------------------------------------------------------
  !> Both direct solvers and FMM calculate the Hartree potential via 
  !! direct additions, without solving the Poisson equation itself.
  !! Direct solvers in any dimension does not require any initialization.
  !! However, fmm requires initialization because, in contrast to Octopus`
  !! direct solvers, FMM reads some parameters from the inp file.
  subroutine poisson_fmm_solve(this, der, pot, rho)  
    type(poisson_fmm_t),         intent(in)    :: this
    type(derivatives_t), target, intent(in)    :: der
    FLOAT,                       intent(inout) :: pot(:)
    FLOAT,                       intent(in)    :: rho(:)

#ifdef HAVE_LIBFM
    integer(8) :: totalcharges

    real(kind = fcs_real_kind_isoc), allocatable :: q(:)  
    real(kind = fcs_real_kind_isoc), allocatable :: pot_lib_fmm(:)
    real(kind = fcs_real_kind_isoc), allocatable :: xyz(:)
    real(kind = fcs_real_kind_isoc), allocatable :: fields(:)
    integer(kind = fcs_integer_kind_isoc)        :: index
    FLOAT,   allocatable :: rho_tmp(:), pot_tmp(:)
    real(8) :: aux
    integer :: ii, jj, ip, gip
    type(c_ptr) ::  ret
    integer, allocatable :: ix(:)
    type(mesh_t), pointer :: mesh

    type(profile_t), save :: poisson_prof, prof_fmm_lib, prof_fmm_corr, prof_fmm_gat

    PUSH_SUB(poisson_fmm_solve)

    mesh => der%mesh

    call profiling_in(poisson_prof, "POISSON_FMM")

    ! allocate buffers.
    SAFE_ALLOCATE(q(this%sp:this%ep))
    SAFE_ALLOCATE(xyz(1:3 * (this%nlocalcharges)))
    SAFE_ALLOCATE(pot_lib_fmm(this%sp:this%ep))
    SAFE_ALLOCATE(fields(1:3 * (this%nlocalcharges)))

    totalcharges = mesh%np_global 

    if (.not. mesh%use_curvilinear) then
      do ii = this%sp, this%ep
        q(ii) = rho(ii) * mesh%vol_pp(1)
      end do
    else
      do ii = this%sp, this%ep
        q(ii) = rho(ii) * mesh%vol_pp(ii)
      end do
    end if

    ! invert the indices
    index = 0
    do ii = this%sp, this%ep
      do jj = 1, der%dim
         index = index + 1 
         xyz(index) = mesh%x(ii, jj)
      end do
    end do

    pot_lib_fmm = M_ZERO 
    fields = M_ZERO

    call profiling_in(prof_fmm_lib, "FMM_LIB")
    ret = fcs_tune(this%handle, this%nlocalcharges, this%nlocalcharges, xyz(1), q(this%sp))
    ret = fcs_run(this%handle, this%nlocalcharges, this%nlocalcharges, xyz(1), q(this%sp), fields(1), &
         pot_lib_fmm(this%sp))
    call profiling_out(prof_fmm_lib)

    call profiling_in(prof_fmm_gat, "FMM_GATHER")
    if (mpi_world%size > 1) then
      !now we need to allgather the results between "states"
#ifdef HAVE_MPI
      call MPI_Allgatherv(pot_lib_fmm(this%sp), this%nlocalcharges, MPI_FLOAT, &
        pot(1), this%dsize(1), this%disps(1), MPI_FLOAT, &
        this%perp_grp%comm, mpi_err)
#endif
    else
      pot = pot_lib_fmm
    end if
    call profiling_out(prof_fmm_gat)

    ! Correction for cell self-interaction in 2D (cell assumed to be circular)
    if (mesh%sb%dim == 2) then
      aux = M_TWO*M_PI * mesh%spacing(1)
      do ii = 1, mesh%np
        pot(ii) = pot(ii) + aux*rho(ii)
      end do
    end if

    SAFE_DEALLOCATE_A(q)
    SAFE_DEALLOCATE_A(xyz)
    SAFE_DEALLOCATE_A(pot_lib_fmm)

    ! Apply the parallel correction
    call profiling_in(prof_fmm_corr, "FMM_CORR")

    SAFE_ALLOCATE(ix(1:mesh%sb%dim))
    SAFE_ALLOCATE(rho_tmp(1:mesh%np_part))
    SAFE_ALLOCATE(pot_tmp(1:mesh%np)) 

    call lalg_copy(mesh%np, rho, rho_tmp)

    ! FMM just calculates contributions from other cells. for
    ! self-interaction cell integration, we include (as traditional in
    ! octopus) an approximate integration using a spherical cell whose
    ! volume is the volume of the actual cell Next line is only valid
    ! for 3D
    if (mesh%sb%dim == 3) then
      if (.not. mesh%use_curvilinear .and. (mesh%spacing(1) == mesh%spacing(2)) .and. &
        (mesh%spacing(2) == mesh%spacing(3)) .and. &
        (mesh%spacing(1) == mesh%spacing(3))) then

        call dderivatives_perform(this%corrector, der, rho_tmp, pot_tmp)

        do ip = 1, mesh%np
          pot(ip) = pot(ip) + pot_tmp(ip)
        end do

      else ! Not common mesh; we add the self-interaction of the cell
        do ii = 1, mesh%np 
          aux = M_TWO*M_PI*(CNST(3.0)*mesh%vol_pp(ii)/(M_PI*CNST(4.0)))**(CNST(2.0)/CNST(3.0))
          pot(ii) = pot(ii) + aux*rho_tmp(ii)
        end do
      end if
    end if

    SAFE_DEALLOCATE_A(ix)
    SAFE_DEALLOCATE_A(rho_tmp)
    SAFE_DEALLOCATE_A(pot_tmp)

    call profiling_out(prof_fmm_corr)
    call profiling_out(poisson_prof)

    POP_SUB(poisson_fmm_solve)
#endif
  end subroutine poisson_fmm_solve

end module poisson_fmm_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
