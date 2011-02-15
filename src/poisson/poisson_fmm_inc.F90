!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
!! $Id: poisson_fmm.F90 7XXX 2011-01-15 22:19:34Z pgarciarisueno and jalberdi$

! ---------------------------------------------------------
subroutine poisson_fmm_init(params_fmm, all_nodes_comm)
  type(poisson_fmm_t), intent(out)   :: params_fmm
  integer,             intent(in)    :: all_nodes_comm

#ifdef HAVE_LIBFM
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
  !% 1 = absolute deltaE error. The error (deltaE) is a fraction of the unity of energy
  !% 2 = relative deltaE error. The error is the given ratio (deltaE) of the total energy
  !% > - Could you explain me what is the difference between considering relative
  !% > or absolute error in the calculations, and why you choose your default as
  !%  > deltaE=E-3, absrel=relative?
  !% The default is just standard error, which fits most situations. It
  !% means, your energy has three significant digits. Lets say the energy of
  !% your system is 1000.0, then the FMM will compute results with a
  !% precision of +-1.
  !% So the result will be
  !% energy=999 ... 1001.
  !% If you change deltaE to 10^-6 it would be something in between
  !% energy=999.999 ... 1000.001
  !% If you do know the magnitude of your energy and set absrel to an
  !% absolute error the situation is different. Setting deltaE to 10^-2 means
  !% you will an energy=999.99..1000.01 which corresponds to 10^⁻5 as
  !% relative error.
  !% Which one you choose is up to you. Since you want to calculate periodic
  !% systems, you may experience very precise results even if you set deltaE
  !% very low. Is is a side effect from the periodicity (totalcharge=0), but
  !% should not bother you at all. You get this kind of extra precision for free.
  !%End
  call parse_integer(datasets_check('AbsRelFMM'), 2, params_fmm%abs_rel_fmm)

  !%Variable DipoleCorrection 
  !%Type integer
  !%Default 0
  !%Section Hamiltonian::Poisson 
  !%Description
  !% Extrinsic/Intrinsic potential
  !%  0 = FMM decides, if correction should be applied
  !%  1 = apply dipole correction  for this run
  !% -1 = disables dipole correction for this run
  !% if you want to compare to classical Ewald use 0 or 1
  !%End
  call parse_integer(datasets_check('DipoleCorrection'), 0, params_fmm%dipole_correction)

  ! FMM: Variable periodic sets periodicity
  ! 0 = open system
  ! 1 = 1D periodic system
  ! 2 = 2D periodic system
  ! 3 = 3D periodic system
  
  call mpi_grp_init(params_fmm%mpi_grp, all_nodes_comm)

  call fmm_init()

  POP_SUB(poisson_fmm_init)
#endif
end subroutine poisson_fmm_init

! ---------------------------------------------------------
subroutine poisson_fmm_end()
#ifdef HAVE_LIBFM
  call fmm_finalize()
#endif
end subroutine poisson_fmm_end

! ---------------------------------------------------------
!> Both direct solvers and FMM calculate the Hartree potential via 
!! direct additions, without solving the Poisson equation itself.
!! Direct solvers in any dimension does not require any initialization.
!! However, fmm requires initialization because, in contrast to Octopus'
!! direct solvers, FMM reads some parameters from the inp file.
subroutine poisson_fmm_solve(this, pot, rho)  
  type(poisson_t),     intent(inout) :: this
  FLOAT,               intent(inout) :: pot(:)
  FLOAT,               intent(in)    :: rho(:)

#ifdef HAVE_LIBFM
  integer(8) :: nlocalcharges
  integer(8) :: totalcharges
  integer(8) :: periodic
  integer(8) :: periodicaxes  !< is always 1
  integer(8) :: dipolecorrection
  integer(8) :: absrel

  real(8), allocatable :: q(:)  
  real(8), allocatable :: potLibFMM(:)
  real(8), allocatable :: xyz(:, :)
  real(8) :: deltaE 
  real(8) :: energyfmm   !< We don't use it, but we cannot remove energyfmm by the moment
  real(8) :: periodlength
  real(8) :: st, en
  real(8) :: aux
  integer :: ii, jj, ierr
  integer :: coords(1:2), sizes(1:2)
  logical :: periods(1:2)
  integer, allocatable :: dstart(:), dend(:)

  PUSH_SUB(poisson_fmm_solve)

  call profiling_in(poisson_prof, "POISSON_FMM")
  call MPI_Cart_get (this%fmm_params%mpi_grp%comm, 2, sizes(1), periods(1), coords(1), mpi_err)

  ASSERT(sizes(1) == this%der%mesh%mpi_grp%size)
  ASSERT(coords(1) == this%der%mesh%mpi_grp%rank)
  
  this%fmm_params%periodic = this%der%mesh%sb%periodic_dim

  if(this%fmm_params%periodic /= 0) then
    if ((this%der%mesh%sb%box_shape == PARALLELEPIPED).and.((this%der%mesh%sb%lsize(1)==this%der%mesh%sb%lsize(2)).and.&
      (this%der%mesh%sb%lsize(1)==this%der%mesh%sb%lsize(3)).and.&
      (this%der%mesh%sb%lsize(2)==this%der%mesh%sb%lsize(3)))) then
      this%fmm_params%periodic_length = this%der%mesh%sb%lsize(1)
    else
      message(1) = "At present, FMM solver for Hartree potential can only deal with cubic boxes. "
      message(2) = " Please, change your Poisson solver or the size or dimensions of your box. "
      call write_fatal(2)
    endif
  endif

  ! set the number of all charges in our system
  nlocalcharges = this%der%mesh%np

  ! allocate buffers.
  SAFE_ALLOCATE(q(nlocalcharges))
  SAFE_ALLOCATE(xyz(1:3,nlocalcharges))
  SAFE_ALLOCATE(potLibFMM(nlocalcharges)) 

  totalcharges = this%der%mesh%np_global 

  !! Beware of this, this is only valid if case(curv_method_uniform); for other cases, vol_pp is to be used
  !! in this case, vol_pp has no meaning, their values are ridiculous
  !! urvilinear coordinates must be included

  if (.not. this%der%mesh%use_curvilinear) then
    aux=this%der%mesh%spacing(1)*this%der%mesh%spacing(2)*this%der%mesh%spacing(3)
    do ii = 1, this%der%mesh%np
      q(ii) = rho(ii)*aux
    end do
  else
    do ii = 1, this%der%mesh%np 
      q(ii) = rho(ii)*this%der%mesh%vol_pp(ii)
    end do
  end if

  ! invert the indexes
  do ii = 1, this%der%mesh%np
    do jj=1, MAX_DIM
      xyz(jj,ii) = this%der%mesh%x(ii,jj)
    end do
  end do

  absrel = this%fmm_params%abs_rel_fmm
  deltaE = this%fmm_params%delta_E_fmm
  potLibFMM = M_ZERO 
  this%fmm_params%periodic=this%der%mesh%sb%periodic_dim
  this%fmm_params%periodic_length=2.0*this%der%mesh%sb%lsize(1)
  periodic = this%fmm_params%periodic
  periodicaxes = 1 
  periodlength = this%fmm_params%periodic_length
  dipolecorrection = this%fmm_params%dipole_correction

  call fmm(totalcharges, nlocalcharges, q, xyz, absrel, deltaE, energyfmm, potLibFMM, &
    periodic, periodicaxes, periodlength, dipolecorrection)

  pot = potLibFMM
  potLibFMM = M_ZERO

  ! FMM just calculates contributions from other cells. for self-interaction cell integration, we include 
  ! (as traditional in octopus) an approximate integration using a spherical cell whose volume is the volume of the actual cell
  ! Next line is only valid for 3D
  if (this%der%mesh%sb%dim==3) then
    if (.not. this%der%mesh%use_curvilinear .and. (this%der%mesh%spacing(1)==this%der%mesh%spacing(2)) .and. &
      (this%der%mesh%spacing(2)==this%der%mesh%spacing(3)) .and. &
      (this%der%mesh%spacing(1)==this%der%mesh%spacing(3))) then
      aux = CNST(2.380077363979553356918)*(this%der%mesh%spacing(1)*this%der%mesh%spacing(2)) 
      do ii = 1, this%der%mesh%np
        pot(ii)=pot(ii)+aux*rho(ii)
      end do
    else
      do ii = 1, this%der%mesh%np 
        aux = M_TWO*M_PI*(3.*this%der%mesh%vol_pp(ii)/(M_PI*4.))**(2./3.)
        pot(ii)=pot(ii)+aux*rho(ii)
      end do
    end if
  end if

  if (this%der%mesh%sb%dim==2) then
    aux = M_TWO*M_PI*this%der%mesh%spacing(1)
    do ii = 1, this%der%mesh%np
      pot(ii)=pot(ii)+aux*rho(ii)
    end do
  end if

  SAFE_DEALLOCATE_A(q)
  SAFE_DEALLOCATE_A(xyz)
  SAFE_DEALLOCATE_A(potLibFMM)

  call profiling_out(poisson_prof)

  POP_SUB(poisson_fmm)
#endif
end subroutine poisson_fmm_solve

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:


