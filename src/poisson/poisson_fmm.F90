!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch, 
!! J. Alberdi, P. Garcia RisueÃ±o, M. Oliveira
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
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
!! 02111-1307, USA.
!!
!! $Id$

#include "global.h"
#ifdef HAVE_LIBFM
#include "fcs_fconfig.h"
#endif

module poisson_fmm_m
  use boundaries_m
  use cube_m
  use datasets_m
  use derivatives_m
  use fft_m
  use geometry_m
  use global_m
  use index_m
  use io_m
  use io_function_m
  use lalg_basic_m
  use loct_math_m
  use math_m
  use mesh_m
  use mesh_cube_parallel_map_m
  use mesh_function_m
  use messages_m
  use mpi_m
  use multicomm_m
  use nl_operator_m
  use par_vec_m
  use parser_m
  use profiling_m
  use simul_box_m
  use stencil_star_m
  use varinfo_m

#ifdef HAVE_LIBFM
  use fcs_module
#else
#define fcs_integer_kind_isoc 4
#endif
  use iso_fortran_env
  use iso_c_binding
    
  implicit none

  private

  public ::                      &
    poisson_fmm_t,               &
    poisson_fmm_init,            &
    poisson_fmm_end,             &
    poisson_fmm_solve

  type poisson_fmm_t
    FLOAT   :: delta_E_fmm
    integer :: abs_rel_fmm
    integer :: dipole_correction
    FLOAT   :: alpha_fmm  !< Alpha for the correction of the FMM
    type(mpi_grp_t) :: all_nodes_grp !< The communicator for all nodes.
    type(mpi_grp_t) :: perp_grp      !< The communicator perpendicular to the mesh communicator.
    integer(kind = fcs_integer_kind_isoc) :: nlocalcharges
    integer    :: sp !< Local start point
    integer    :: ep !< Local end point
    integer, pointer :: disps(:)
    integer, pointer :: dsize(:) !< Local size
    type(nl_operator_t) :: corrector
    type(derivatives_t), pointer :: der
    type(c_ptr) ::  handle !< The FMM identifier
  end type poisson_fmm_t

contains

  ! ---------------------------------------------------------
  !> Initialises the FMM parameters and vectors. Also it calls to 
  !! the library initialisation.
  subroutine poisson_fmm_init(this, der, all_nodes_comm)
    type(poisson_fmm_t),         intent(out)   :: this
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

    real(kind = fcs_real_kind), parameter           ::  BOX_SIZE = 1.00d0
    logical ::  short_range_flag = .true.
    real(kind = fcs_real_kind_isoc), dimension(3)   ::  box_a !! = (/BOX_SIZE,0.0d0,0.0d0/)
    real(kind = fcs_real_kind_isoc), dimension(3)   ::  box_b !!= (/0.0d0,BOX_SIZE,0.0d0/)
    real(kind = fcs_real_kind_isoc), dimension(3)   ::  box_c !! = (/0.0d0,0.0d0,BOX_SIZE/)
    real(kind = fcs_real_kind_isoc), dimension(3)   ::  offset = (/M_ZERO,M_ZERO,M_ZERO/)
    logical, dimension(3)                           ::  periodicity = (/.false.,.false.,.false./)
    integer(kind = fcs_integer_kind_isoc)           ::  total_particles

    integer(kind = fcs_integer_kind_isoc)           ::  local_particle_count = -1
    real(kind = fcs_real_kind_isoc), dimension(8)   ::  local_charges
    real(kind = fcs_real_kind_isoc), dimension(24)  ::  local_coordinates

    integer(8) :: periodic
    integer(8) :: periodicaxes  !< is always 1
    real(8) :: periodic_length

    PUSH_SUB(poisson_fmm_init)

    method = "fmm"
    !%Variable DeltaEFMM
    !%Type float
    !%Default 0.0001 
    !%Section Hamiltonian::Poisson 
    !%Description
    !% Parameter for absolute or relative convergence of FMM.
    !% Sets energy error bound.
    !% Strong inhomogeneous systems may violate the error bound.
    !% For inhomogeneous systems we have an error-controlled sequential version available
    !% (from Ivo Kabadshow).
    !%End
    call parse_float(datasets_check('DeltaEFMM'), CNST(1e-4), this%delta_E_fmm)

    !%Variable AbsRelFMM 
    !%Type integer
    !%Default 2
    !%Section Hamiltonian::Poisson 
    !%Description
    !% Sets type of error bound.
    !% 0 = 10^-3 relative error.
    !% 1 = absolute delta_E error. The error (delta_E) is a fraction of the unity of energy
    !% 2 = relative delta_E error. The error is the given ratio (delta_E) of the total energy
    !% > - Could you explain me what is the difference between considering relative
    !% > or absolute error in the calculations, and why you choose your default as
    !%  > delta_E=E-3, absrel=relative?
    !% The default is just standard error, which fits most situations. It
    !% means, your energy has three significant digits. Lets say the energy of
    !% your system is 1000.0, then the FMM will compute results with a
    !% precision of +-1.
    !% So the result will be
    !% energy=999 ... 1001.
    !% If you change delta_E to 10^-6 it would be something in between
    !% energy=999.999 ... 1000.001
    !% If you do know the magnitude of your energy and set absrel to an
    !% absolute error the situation is different. Setting delta_E to 10^-2 means
    !% you will an energy=999.99..1000.01 which corresponds to 10^5 as
    !% relative error.
    !% Which one you choose is up to you. Since you want to calculate periodic
    !% systems, you may experience very precise results even if you set delta_E
    !% very low. It is a side effect from the periodicity (totalcharge=0), but
    !% should not bother you at all. You get this kind of extra precision for free.
    !%End
    call parse_integer(datasets_check('AbsRelFMM'), 2, this%abs_rel_fmm)

    !%Variable DipoleCorrection 
    !%Type integer
    !%Default 0
    !%Section Hamiltonian::Poisson 
    !%Description
    !% Extrinsic/Intrinsic potential.
    !% If you want to compare to classical Ewald use 0 or 1.
    !%Option 0
    !% FMM decides whether correction should be applied.
    !%Option 1
    !% Apply dipole correction.
    !%Option -1
    !% Disables dipole correction.
    !%End
    call parse_integer(datasets_check('DipoleCorrection'), 0, this%dipole_correction)

    !%Variable AlphaFMM
    !%Type float
    !%Default 0.291262136
    !%Section Hamiltonian::Poisson 
    !%Description
    !% Parameter for the correction of the self-interaction of the
    !% electrostatic Hartree potential. The default value is 0.291262136.
    !%End
    call parse_float(datasets_check('AlphaFMM'), CNST(0.291262136), this%alpha_fmm)

    ! FMM: Variable periodic sets periodicity
    ! 0 = open system
    ! 1 = 1D periodic system
    ! 2 = 2D periodic system
    ! 3 = 3D periodic system

    call mpi_grp_init(this%all_nodes_grp, all_nodes_comm)

    this%der => der

    if (mpi_world%size == 1) then
      cdim = 1

      SAFE_ALLOCATE(this%disps(1))
      SAFE_ALLOCATE(dend(1))
      SAFE_ALLOCATE(this%dsize(1))

      dend = der%mesh%np
      this%sp = 1 
      this%ep = der%mesh%np
      this%dsize(1) = der%mesh%np
      this%disps = 0
      this%nlocalcharges = this%dsize(1)

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
    periodic = mesh%sb%periodic_dim

    if(periodic /= 0) then
      if ((mesh%sb%box_shape == PARALLELEPIPED) .and. ((mesh%sb%lsize(1) == mesh%sb%lsize(2)) .and. &
        (mesh%sb%lsize(1) == mesh%sb%lsize(3)) .and. &
        (mesh%sb%lsize(2) == mesh%sb%lsize(3)))) then
        periodic_length = mesh%sb%lsize(1)
      else
        message(1) = "At present, FMM solver for Hartree potential can only deal with cubic boxes. "
        message(2) = " Please, change your Poisson solver or the size or dimensions of your box. "
        call messages_fatal(2)
      end if
    end if

    total_particles = mesh%np_part_global !! mesh%np_globla
    total_particles = mesh%np_global
    ret = fcs_init(this%handle, trim(adjustl(method)) // c_null_char, subcomm)

    box_a = (/mesh%sb%lsize(1)*2+1,M_ZERO,M_ZERO/)
    box_b = (/M_ZERO,mesh%sb%lsize(2)*2+1,M_ZERO/)
    box_c = (/M_ZERO,M_ZERO,mesh%sb%lsize(3)*2+1/)

    ret = fcs_common_set(this%handle, short_range_flag, box_a, box_b, box_c, offset, periodicity, total_particles)

    ! this is how you set a relative error in scafacos for the FMM
    ret = fcs_fmm_set_tolerance_energy( this%handle, this%delta_E_fmm );

    call nl_operator_init(this%corrector, "FMM Correction")
    call stencil_star_get_lapl(this%corrector%stencil, der%mesh%sb%dim, 2)
    call nl_operator_build(this%der%mesh, this%corrector, der%mesh%np, const_w = .not. this%der%mesh%use_curvilinear)

    do is = 1, this%corrector%stencil%size

      select case(sum(abs(this%corrector%stencil%points(1:MAX_DIM, is))))
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

    call nl_operator_update_weights(this%corrector)
    POP_SUB(poisson_fmm_init)
#endif
  end subroutine poisson_fmm_init

  ! ---------------------------------------------------------
  !> Release memory and call to end the library
  subroutine poisson_fmm_end(this)
    type(poisson_fmm_t), intent(inout) :: this
    type(c_ptr) ::  ret
    
#ifdef HAVE_LIBFM
    PUSH_SUB(poisson_fmm_end)

    call nl_operator_end(this%corrector)

    if (mpi_world%size > 1) call MPI_Comm_free(this%perp_grp%comm, mpi_err)
    SAFE_DEALLOCATE_P(this%disps)
    SAFE_DEALLOCATE_P(this%dsize)
!!$    ret = fcs_destroy(this%handle)

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
    type(poisson_fmm_t),         intent(inout) :: this
    type(derivatives_t), target, intent(in)    :: der
    FLOAT,                       intent(inout) :: pot(:)
    FLOAT,                       intent(in)    :: rho(:)

#ifdef HAVE_LIBFM
    integer(8) :: totalcharges
    integer(8) :: periodic
    integer(8) :: periodicaxes  !< is always 1
    integer(8) :: dipolecorrection
    integer(8) :: abs_rel

    real(kind = fcs_real_kind_isoc), allocatable :: q(:)  
    real(kind = fcs_real_kind_isoc), allocatable :: pot_lib_fmm(:)
    real(kind = fcs_real_kind_isoc), allocatable :: xyz(:)
    FLOAT,   allocatable :: rho_tmp(:), pot_tmp(:)
    real(8) :: delta_E 
    real(8) :: energy_fmm   !< We don`t use it, but we cannot remove energy_fmm for the moment
    real(8) :: periodic_length
    real(8) :: aux
    integer :: ii, jj, ip, gip
    type(c_ptr) ::  ret
    real(kind = fcs_real_kind_isoc), allocatable         ::  fields(:)
    integer(kind = fcs_integer_kind_isoc)           ::  index
    integer, allocatable :: ix(:)
    FLOAT :: aux1
    type(mesh_t), pointer :: mesh

    type(profile_t), save :: poisson_prof, prof_fmm_lib, prof_fmm_corr, prof_fmm_gat

    PUSH_SUB(poisson_fmm_solve)

    mesh => der%mesh

    call profiling_in(poisson_prof, "POISSON_FMM")

    ! allocate buffers.
    SAFE_ALLOCATE(q(this%sp:this%ep))
    SAFE_ALLOCATE(xyz(3 * (this%nlocalcharges)))
    SAFE_ALLOCATE(pot_lib_fmm(this%sp:this%ep))
    SAFE_ALLOCATE(fields(3 * (this%nlocalcharges)))

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
      do jj = 1, MAX_DIM
         index = index + 1 
         xyz(index) = mesh%x(ii, jj)
      end do
    end do

    abs_rel = this%abs_rel_fmm
    delta_E = this%delta_E_fmm
    pot_lib_fmm = M_ZERO 
    fields = M_ZERO
    periodic_length = CNST(2.0)*mesh%sb%lsize(1)
    periodicaxes = 1
    dipolecorrection = this%dipole_correction

    call profiling_in(prof_fmm_lib, "FMM_LIB")

    ret = fcs_tune(this%handle, this%nlocalcharges, this%nlocalcharges, xyz(1), q(this%sp))

    ret = fcs_run(this%handle, this%nlocalcharges, this%nlocalcharges, xyz(1), q(this%sp), fields(1), &
         pot_lib_fmm(this%sp))

    call profiling_out(prof_fmm_lib)

    call profiling_in(prof_fmm_gat, "FMM_GATHER")
    if (mpi_world%size > 1) then
      !now we need to allgather the results between "states"
      call MPI_Allgatherv(pot_lib_fmm(this%sp), this%nlocalcharges, MPI_FLOAT, &
        pot(1), this%dsize(1), this%disps(1), MPI_FLOAT, &
        this%perp_grp%comm, mpi_err)
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

    ! FMM just calculates contributions from other cells. for self-interaction cell integration, we include 
    ! (as traditional in octopus) an approximate integration using a spherical cell whose volume is the volume of the actual cell
    ! Next line is only valid for 3D
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

end module poisson_fmm_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:

