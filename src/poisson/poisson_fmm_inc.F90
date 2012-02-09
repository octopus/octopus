!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch, 
!! J. Alberdi, P. Garcia Risueño, M. Oliveira
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

! ---------------------------------------------------------
!> Initialises the FMM parameters and vectors. Also it calls to 
!! the library initialisation.
subroutine poisson_fmm_init(params_fmm, mesh, all_nodes_comm)
  type(poisson_fmm_t), intent(out)   :: params_fmm
  type(mesh_t),        intent(in)    :: mesh
  integer,             intent(in)    :: all_nodes_comm

#ifdef HAVE_LIBFM
  logical, allocatable :: remains(:)
  integer, allocatable :: dend(:)
  integer :: subcomm, cdim

  PUSH_SUB(poisson_fmm_init)

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
  call parse_float(datasets_check('DeltaEFMM'), CNST(1e-4), params_fmm%delta_E_fmm)

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
  call parse_integer(datasets_check('AbsRelFMM'), 2, params_fmm%abs_rel_fmm)

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
  call parse_integer(datasets_check('DipoleCorrection'), 0, params_fmm%dipole_correction)
  
  !%Variable AlphaFMM
  !%Type float
  !%Default 0.291262136
  !%Section Hamiltonian::Poisson 
  !%Description
  !% Parameter for the correction of the self-interaction of the
  !% electrostatic Hartree potential. The default value is 0.291262136.
  !%End
  call parse_float(datasets_check('AlphaFMM'), CNST(0.291262136), params_fmm%alpha_fmm)

  ! FMM: Variable periodic sets periodicity
  ! 0 = open system
  ! 1 = 1D periodic system
  ! 2 = 2D periodic system
  ! 3 = 3D periodic system
  
  call mpi_grp_init(params_fmm%all_nodes_grp, all_nodes_comm)

  if (mpi_world%size == 1) then
    cdim = 1

    SAFE_ALLOCATE(params_fmm%disps(1))
    SAFE_ALLOCATE(dend(1))
    SAFE_ALLOCATE(params_fmm%dsize(1))
    
    dend = mesh%np
    params_fmm%sp = 1 
    params_fmm%ep = mesh%np
    params_fmm%dsize(1) = mesh%np
    params_fmm%disps = 0
    params_fmm%nlocalcharges = params_fmm%dsize(1)
    
  else 

    call MPI_Cartdim_get(params_fmm%all_nodes_grp%comm, cdim, mpi_err)
 
    SAFE_ALLOCATE(remains(1:cdim))
    
    remains = .true.
    remains(1) = .false.
    
    call MPI_Cart_sub(params_fmm%all_nodes_grp%comm, remains(1), subcomm, mpi_err)
    call mpi_grp_init(params_fmm%perp_grp, subcomm)

    SAFE_ALLOCATE(params_fmm%disps(1:params_fmm%perp_grp%size))
    SAFE_ALLOCATE(dend(1:params_fmm%perp_grp%size))
    SAFE_ALLOCATE(params_fmm%dsize(1:params_fmm%perp_grp%size))
    
    call multicomm_divide_range(mesh%np, params_fmm%perp_grp%size, params_fmm%disps, dend, params_fmm%dsize)
    
    params_fmm%sp = params_fmm%disps(params_fmm%perp_grp%rank + 1)
    params_fmm%ep = dend(params_fmm%perp_grp%rank + 1)
    params_fmm%nlocalcharges = params_fmm%dsize(params_fmm%perp_grp%rank + 1)
    params_fmm%disps = params_fmm%disps - 1

  end if
  call fmm_init()

  POP_SUB(poisson_fmm_init)
#endif
end subroutine poisson_fmm_init

! ---------------------------------------------------------
!> Release memory and call to end the library
subroutine poisson_fmm_end(params_fmm)
  type(poisson_fmm_t), intent(inout) :: params_fmm

#ifdef HAVE_LIBFM
  PUSH_SUB(poisson_fmm_end)

  if (mpi_world%size > 1) call MPI_Comm_free(params_fmm%perp_grp%comm, mpi_err)
  SAFE_DEALLOCATE_P(params_fmm%disps)
  SAFE_DEALLOCATE_P(params_fmm%dsize)
  call fmm_finalize()

  POP_SUB(poisson_fmm_end)
#endif
end subroutine poisson_fmm_end

! ---------------------------------------------------------
!> Both direct solvers and FMM calculate the Hartree potential via 
!! direct additions, without solving the Poisson equation itself.
!! Direct solvers in any dimension does not require any initialization.
!! However, fmm requires initialization because, in contrast to Octopus`
!! direct solvers, FMM reads some parameters from the inp file.
subroutine poisson_fmm_solve(this, pot, rho)  
  type(poisson_t),     intent(inout) :: this
  FLOAT,               intent(inout) :: pot(:)
  FLOAT,               intent(in)    :: rho(:)

#ifdef HAVE_LIBFM
  integer(8) :: totalcharges
  integer(8) :: periodic
  integer(8) :: periodicaxes  !< is always 1
  integer(8) :: dipolecorrection
  integer(8) :: abs_rel

  real(8), allocatable :: q(:)  
  real(8), allocatable :: pot_lib_fmm(:)
  real(8), allocatable :: xyz(:, :)
  FLOAT,   allocatable :: rho_tmp(:)
  real(8) :: delta_E 
  real(8) :: energy_fmm   !< We don`t use it, but we cannot remove energy_fmm for the moment
  real(8) :: periodic_length
  real(8) :: aux
  integer :: ii, jj, sp, ep, ip, gip

  integer, allocatable :: ix(:)
  FLOAT :: aux1
  type(mesh_t), pointer :: mesh

  type(profile_t), save :: prof_fmm_lib, prof_fmm_corr, prof_fmm_gat

  PUSH_SUB(poisson_fmm_solve)

  mesh => this%der%mesh
  sp = this%params_fmm%sp
  ep = this%params_fmm%ep
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
    endif
  endif

  call profiling_in(poisson_prof, "POISSON_FMM")

  ! allocate buffers.
  SAFE_ALLOCATE(q(sp:ep))
  SAFE_ALLOCATE(xyz(1:3, sp:ep))
  SAFE_ALLOCATE(pot_lib_fmm(sp:ep)) 

  totalcharges = mesh%np_global 

  if (.not. mesh%use_curvilinear) then
    do ii = sp, ep
      q(ii) = rho(ii) * mesh%vol_pp(1)
    end do
  else
    do ii = sp, ep
      q(ii) = rho(ii) * mesh%vol_pp(ii)
    end do
  end if

  ! invert the indices
  do ii = sp, ep
    do jj = 1, MAX_DIM
      xyz(jj, ii) = mesh%x(ii, jj)
    end do
  end do

  abs_rel = this%params_fmm%abs_rel_fmm
  delta_E = this%params_fmm%delta_E_fmm
  pot_lib_fmm = M_ZERO 
  periodic_length = CNST(2.0)*mesh%sb%lsize(1)
  periodicaxes = 1
  dipolecorrection = this%params_fmm%dipole_correction

  call profiling_in(prof_fmm_lib, "FMM_LIB")
  call fmm(totalcharges, this%params_fmm%nlocalcharges, q(sp), xyz(1, sp), abs_rel, delta_E, energy_fmm, &
    pot_lib_fmm(sp), periodic, periodicaxes, periodic_length, dipolecorrection)
  call profiling_out(prof_fmm_lib)
  
  call profiling_in(prof_fmm_gat, "FMM_GATHER")
  if (mpi_world%size > 1) then
    !now we need to allgather the results between "states"
    call MPI_Allgatherv(pot_lib_fmm(sp), this%params_fmm%nlocalcharges, MPI_FLOAT, &
         pot(1), this%params_fmm%dsize(1), this%params_fmm%disps(1), MPI_FLOAT, &
         this%params_fmm%perp_grp%comm, mpi_err)
  else
    pot = pot_lib_fmm
  end if
  call profiling_out(prof_fmm_gat)

  ! Correction for cell self-interaction in 2D (cell assumed to be circular)
  if (mesh%sb%dim == 2) then
    aux = M_TWO*M_PI * mesh%spacing(1)
    do ii = 1, mesh%np
      pot(ii) = pot(ii) + aux * rho(ii)
    end do
  end if

  SAFE_DEALLOCATE_A(q)
  SAFE_DEALLOCATE_A(xyz)
  SAFE_DEALLOCATE_A(pot_lib_fmm)

  ! Apply the parallel correction
  call profiling_in(prof_fmm_corr, "FMM_CORR")
  SAFE_ALLOCATE(ix(1:mesh%sb%dim))
  SAFE_ALLOCATE(rho_tmp(1:mesh%np_part))
  do ip = 1, mesh%np
    rho_tmp(ip) = rho(ip)
  end do
  do ip = mesh%np + 1, mesh%np_part
    rho_tmp(ip) = M_ZERO
  end do
  if (mesh%parallel_in_domains) then
#ifdef HAVE_MPI
    call dvec_ghost_update(mesh%vp, rho_tmp)
#endif
  end if
  
  ! FMM just calculates contributions from other cells. for self-interaction cell integration, we include 
  ! (as traditional in octopus) an approximate integration using a spherical cell whose volume is the volume of the actual cell
  ! Next line is only valid for 3D
  if (mesh%sb%dim == 3) then
    if (.not. mesh%use_curvilinear .and. (mesh%spacing(1) == mesh%spacing(2)) .and. &
      (mesh%spacing(2) == mesh%spacing(3)) .and. &
      (mesh%spacing(1) == mesh%spacing(3))) then

      ! Corrections for first neighbours obtained with linear interpolation      
      ! First we obtain the densities in neighbouring points ip+1/2
      ! Iterate over all local points of rho (1 to mesh%np)
      do ip = 1, mesh%np
        if (mesh%parallel_in_domains) then
#ifdef HAVE_MPI
          ! Get the global point from the local point
          gip = mesh%vp%local(mesh%vp%xlocal(mesh%vp%partno) + ip - 1)
#endif
        else
          gip = ip
        end if

        ! Get x,y,z indices of the global point
        call index_to_coords(mesh%idx, mesh%sb%dim, gip, ix)

        !!Correction for FMM with semi-neighbours (terms are merged for computational efficiency)

        aux1 = M_ZERO  

        aux1 = aux1 - rho_tmp(vec_index2local(mesh%vp,mesh%idx, ix, 1, -2)) &
             - rho_tmp(vec_index2local(mesh%vp,mesh%idx, ix, 1, 2)) &
             - rho_tmp(vec_index2local(mesh%vp,mesh%idx, ix, 2, -2)) - rho_tmp(vec_index2local(mesh%vp,mesh%idx, ix, 2, 2)) &
             - rho_tmp(vec_index2local(mesh%vp,mesh%idx, ix, 3, -2)) - rho_tmp(vec_index2local(mesh%vp,mesh%idx, ix, 3, 2)) 

        aux1 = aux1/M_FOUR

        aux1 = aux1 + rho_tmp(vec_index2local(mesh%vp,mesh%idx, ix, 1, -1)) &
             + rho_tmp(vec_index2local(mesh%vp,mesh%idx, ix, 1, 1)) &
             + rho_tmp(vec_index2local(mesh%vp,mesh%idx, ix, 2, -1)) + rho_tmp(vec_index2local(mesh%vp,mesh%idx, ix, 2, 1)) &
             + rho_tmp(vec_index2local(mesh%vp,mesh%idx, ix, 3, -1)) + rho_tmp(vec_index2local(mesh%vp,mesh%idx, ix, 3, 1)) 

        aux1 = aux1/16.0
        
        aux1 = aux1 + rho_tmp(ip) * (27.0/32.0 + (M_ONE - this%params_fmm%alpha_fmm) * M_TWO * M_PI * (3./(M_PI*4.))**(2./3.))

        aux1 = aux1 * (mesh%spacing(1) * mesh%spacing(2))

        ! Apply the correction to the potential
        pot(ip) = pot(ip) + aux1

      end do

    else ! Not common mesh; we add the self-interaction of the cell
      do ii = 1, mesh%np 
        aux = M_TWO * M_PI * (3. * mesh%vol_pp(ii)/(M_PI * 4.))**(2./3.)
        pot(ii) = pot(ii) + aux * rho_tmp(ii)
      end do
    end if
  end if

  SAFE_DEALLOCATE_A(ix)
  SAFE_DEALLOCATE_A(rho_tmp)

  call profiling_out(prof_fmm_corr)
  call profiling_out(poisson_prof)
  
  POP_SUB(poisson_fmm_solve)
#endif
end subroutine poisson_fmm_solve

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
