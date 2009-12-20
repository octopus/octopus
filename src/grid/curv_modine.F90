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
!! $Id$

#include "global.h"

! This module implements the curvilinear coordinates given in
! N. A. Modine, G. Zumbach, and E. Kaxiras, Phys. Rev. B 55, 10289-10301 (1997) 
!
! The local refinement was changed for a simple exponential.
! I believe that the recipe given by the authors is too complicated
! for me to sort it out.

module curv_modine_m
  use datasets_m
  use geometry_m
  use geometry_m
  use global_m
  use parser_m
  use messages_m
  use root_solver_m
  use profiling_m
  use simul_box_m
  use unit_m
  use unit_system_m

  implicit none

  private
  public ::                  &
    curv_modine_t,           &
    curv_modine_init,        &
    curv_modine_end,         &
    curv_modine_chi2x,       &
    curv_modine_jacobian_inv

  type curv_modine_t
    FLOAT :: L(MAX_DIM)    ! size of the box
    FLOAT :: xbar          ! size of central flat region (in units of L)
    FLOAT :: Jbar          ! increase in density of points is 1/J

    FLOAT,            pointer :: Jlocal(:)  ! local (around the atoms) refinement
    FLOAT,            pointer :: Jrange(:)  ! local refinement range

    FLOAT, pointer :: chi_atoms(:,:)
    FLOAT, pointer :: csi(:,:)

    integer :: natoms
  end type curv_modine_t

  integer, parameter :: qq = 3

  type(simul_box_t),   pointer :: sb_p
  type(curv_modine_t), pointer :: cv_p
  FLOAT,           allocatable :: x_p(:)

contains

  ! ---------------------------------------------------------
  subroutine getf2(csi, f, jf)
    FLOAT, intent(in)  :: csi(:)
    FLOAT, intent(out) :: f(:), jf(:, :)

    integer :: i1, j1, i2, j2, index1, index2
    FLOAT :: x(MAX_DIM), chi2(MAX_DIM), rr, dd, dd2

    ! first we fill in cv%csi with the values we have
    index1 = 1
    do i1 = 1, cv_p%natoms
      do j1 = 1, sb_p%dim
        cv_p%csi(j1, i1) = csi(index1)
        index1 = index1 + 1
      end do
    end do

    ! get f and jf
    jf(:,:) = M_ZERO
    do i1 = 1, cv_p%natoms
      call curv_modine_chi2chi2(sb_p, cv_p, cv_p%chi_atoms(:,i1), chi2)
      x(:) = chi2(:)

      do i2 = 1, cv_p%natoms
        rr = sqrt(sum((chi2(:) - cv_p%csi(:,i2))**2))
        dd = exp(-rr**2/(M_TWO*cv_p%Jrange(i2)**2))

        x(:) = x(:) - cv_p%Jlocal(i2)*(chi2(:) - cv_p%csi(:,i2)) * dd
      end do

      do j1 = 1, sb_p%dim
        index1 = (i1-1)*sb_p%dim + j1
        f(index1) = x(j1) - x_p(index1)

        do i2 = 1, cv_p%natoms
          rr  = sqrt(sum((chi2 - cv_p%csi(:,i2))**2))
          dd  = exp(-rr**2/(M_TWO*cv_p%Jrange(i2)**2))
          dd2 = -M_TWO/(M_TWO*cv_p%Jrange(i2)**2)*dd

          index2 = (i2-1)*sb_p%dim + j1
          jf(index1, index2) = cv_p%Jlocal(i2) * dd

          do j2 = 1, sb_p%dim
            index2 = (i2-1)*sb_p%dim + j2

            jf(index1, index2) =  jf(index1, index2) + cv_p%Jlocal(i2) * dd2 * &
              (chi2(j1) - cv_p%csi(j1,i2))*(chi2(j2) - cv_p%csi(j2,i2))
          end do
        end do
      end do
    end do

  end subroutine getf2

  ! ---------------------------------------------------------
  subroutine curv_modine_init(sb, geo, cv)
    type(simul_box_t),   target, intent(in)   :: sb
    type(geometry_t),    target, intent(in)   :: geo
    type(curv_modine_t), target, intent(out)  :: cv

    call parse_float(datasets_check('CurvModineXBar'), M_ONE/M_THREE, cv%xbar)
    call parse_float(datasets_check('CurvModineJBar'), M_HALF, cv%Jbar)

    cv%L = M_ZERO
    cv%L(1:sb%dim) = sb%lsize(1:sb%dim) / cv%Jbar

    if(cv%xbar<M_ZERO.or.cv%xbar>M_ONE) then
      message(1) = 'The parameter "CurvModineXBar" must lie between 0 and 1.'
      call write_fatal(1)
    end if

    SAFE_ALLOCATE(cv%Jlocal(1:geo%natoms))
    SAFE_ALLOCATE(cv%Jrange(1:geo%natoms))

    ! WARNING: the reading has to be done for each atom kind
    call parse_float(datasets_check('CurvModineJlocal'), CNST(0.25), cv%Jlocal(1))
    call parse_float(datasets_check('CurvModineJrange'), units_from_atomic(units_inp%length, M_TWO), cv%Jrange(1))

    cv%Jrange(1) = units_to_atomic(units_inp%length, cv%Jrange(1))

    if(cv%Jlocal(1)<M_ZERO.or.cv%Jlocal(1)>M_ONE) then
      message(1) = 'The parameter "CurvModineJlocal" must lie between 0 and 1.'
      call write_fatal(1)
    end if

    cv%Jlocal(:) = cv%Jlocal(1)
    cv%Jrange(:) = cv%Jrange(1)

    call find_atom_points()
    call optimize()

    cv%natoms = geo%natoms

  contains
    subroutine find_atom_points()
      integer :: i, jj

      ! Initialize csi
      SAFE_ALLOCATE(cv%csi(1:sb%dim, 1:geo%natoms))
      do i = 1, geo%natoms
        cv%csi(1:sb%dim,i) = geo%atom(i)%x(1:sb%dim)
      end do

      ! get first estimate for chi_atoms
      SAFE_ALLOCATE(cv%chi_atoms(1:sb%dim, 1:geo%natoms))
      do jj = 1, 10  ! WARNING: make something better
        do i = 1, geo%natoms
          call curv_modine_x2chi(sb, cv, geo%atom(i)%x, cv%chi_atoms(:,i))
        end do
        cv%csi(:,:) = cv%chi_atoms(:,:)
      end do

      do i = 1, geo%natoms
        ! these are the chi positions where we want the atoms
        ! FIXME: sb should not know about the spacing
        cv%chi_atoms(:,i) = nint(cv%chi_atoms(:,i)/sb%spacing(:))*sb%spacing(:)
      end do

    end subroutine find_atom_points

    subroutine optimize()
      logical :: conv
      type(root_solver_t) :: rs
      integer :: i, j, index
      FLOAT, allocatable :: my_csi(:), start_csi(:)

      call root_solver_init(rs, sb%dim, &
        solver_type = ROOT_NEWTON, maxiter = 500, abs_tolerance = CNST(1.0e-10))

      sb_p  => sb
      cv_p  => cv

      SAFE_ALLOCATE(x_p(1:sb%dim*geo%natoms))
      SAFE_ALLOCATE(my_csi(1:sb%dim*geo%natoms))
      SAFE_ALLOCATE(start_csi(1:sb%dim*geo%natoms))

      do i = 1, geo%natoms
        do j = 1, sb%dim
          index = (i-1)*sb%dim + j
          x_p(index)       = geo%atom(i)%x(j)
          start_csi(index) = cv%chi_atoms(j, i)
        end do
      end do

      call droot_solver_run(rs, getf2, my_csi, conv, startval=start_csi)

      if(.not.conv) then
        message(1) = "During the construction of the adaptive grid, the Newton-Raphson"
        message(2) = "method did not converge"
        call write_fatal(2)
      end if

      ! Now set csi to the new values
      do i = 1, geo%natoms
        do j = 1, sb%dim
          index = (i-1)*sb_p%dim + j
          cv_p%csi(j, i) = my_csi(index)
        end do
      end do

      SAFE_DEALLOCATE_A(x_p)
      SAFE_DEALLOCATE_A(my_csi)
      SAFE_DEALLOCATE_A(start_csi)

      nullify(sb_p); nullify(cv_p)

    end subroutine optimize

  end subroutine curv_modine_init


  ! ---------------------------------------------------------
  subroutine curv_modine_end(cv)
    type(curv_modine_t), intent(inout) :: cv

    SAFE_DEALLOCATE_P(cv%Jlocal)
    SAFE_DEALLOCATE_P(cv%Jrange)
    SAFE_DEALLOCATE_P(cv%chi_atoms)

  end subroutine curv_modine_end


  ! ---------------------------------------------------------
  subroutine curv_modine_chi2chi2(sb, cv, chi_, chi2, Jac)
    type(simul_box_t),   intent(in)  :: sb
    type(curv_modine_t), intent(in)  :: cv
    FLOAT,               intent(in)  :: chi_(:)  ! chi_(sb%dim)
    FLOAT,               intent(out) :: chi2(:)  ! chi2(sb%dim)
    FLOAT,     optional, intent(out) :: Jac(:)   ! the Jacobian of this transformation is diagonal

    FLOAT :: chibar(MAX_DIM), rr, chi
    logical :: neg
    integer :: i

    chibar(1:sb%dim) = cv%xbar*cv%L(1:sb%dim)

    do i = 1, sb%dim
      neg = (chi_(i) < 0)
      chi = abs(chi_(i))

      chi2(i)  = cv%Jbar * chi
      if(present(Jac)) Jac(i) = cv%Jbar

      if(chi > chibar(i)) then
        rr = (chi-chibar(i))/(cv%L(i)-chibar(i))

        chi2(i)  = chi2(i) + cv%L(i)/M_TWO*(1-cv%Jbar) * rr**qq *   &
           (qq + M_ONE - (qq - M_ONE)*rr)

        if(present(Jac)) then
          Jac(i) = Jac(i) + cv%L(i)/M_TWO*(1-cv%Jbar) * rr**(qq-1)/(cv%L(i)-chibar(i)) *   &
            (qq*(qq+1) - (qq**2-1)*rr)
        end if
      end if

      if(neg) chi2(i) = -chi2(i)
      ! CHECK if Jacobian does not have to be negated!
    end do

  end subroutine curv_modine_chi2chi2

  ! ---------------------------------------------------------
  subroutine curv_modine_chi2x(sb, cv, chi_, x)
    type(simul_box_t),   intent(in)  :: sb
    type(curv_modine_t), intent(in)  :: cv
    FLOAT,               intent(in)  :: chi_(:)  ! chi_(sb%dim)
    FLOAT,               intent(out) :: x(:)     !   x (sb%dim)

    FLOAT :: chi2(MAX_DIM), rr, dd
    integer :: i

    call curv_modine_chi2chi2(sb, cv, chi_, chi2)

    x(:) = chi2(:)
    do i = 1, cv%natoms
      rr = max(sqrt(sum((chi2(:) - cv%csi(:,i))**2)), CNST(1e-6))
      dd = exp(-rr**2/(M_TWO*cv%Jrange(i)**2))

      x(:) = x(:) - cv%Jlocal(i)*(chi2(:) - cv%csi(:,i)) * dd
    end do

  end subroutine curv_modine_chi2x


  ! ---------------------------------------------------------
  subroutine curv_modine_jacobian_inv(sb, cv, chi_, x, J)
    type(simul_box_t),   intent(in)  :: sb
    type(curv_modine_t), intent(in)  :: cv
    FLOAT,               intent(in)  :: chi_(:)  ! chi(sb%dim)
    FLOAT,               intent(out) :: x(:)     ! x(sb%dim)
    FLOAT,               intent(out) :: J(:,:)   ! J(sb%dim,sb%dim), the Jacobian

    FLOAT :: chi2(MAX_DIM), rr, dd, J2(MAX_DIM)
    integer :: i, ix, iy

    call curv_modine_chi2chi2(sb, cv, chi_, chi2, J2)

    ! initialize both x and the jacobian
    x(:)   = chi2(:)
    J(:,:) = M_ZERO
    do ix = 1, sb%dim
      J(ix, ix) = M_ONE
    end do

    do i = 1, cv%natoms
      rr = max(sqrt(sum(chi2(:) - cv%csi(:,i))**2), CNST(1e-6))
      dd = exp(-rr**2/(M_TWO*cv%Jrange(i)**2))

      x(:) = x(:) -  cv%Jlocal(i)*(chi2(:) - cv%csi(:,i)) * dd

      do ix = 1, sb%dim
        J(ix, ix) = J(ix, ix) - cv%Jlocal(i) * dd
        do iy = 1, sb%dim
          J(ix, iy) = J(ix, iy) + cv%Jlocal(i)*(chi2(ix)-cv%csi(ix,i))*(chi2(iy)-cv%csi(iy,i)) * &
             M_TWO/(M_TWO*cv%Jrange(i)**2) * dd
        end do
      end do
    end do

    do ix = 1, sb%dim
      J(ix,:) = J(ix,:)*J2(:)
    end do

  end subroutine curv_modine_jacobian_inv


  ! ---------------------------------------------------------
  subroutine getf(y, f, jf)
    FLOAT, intent(in)  :: y(:)
    FLOAT, intent(out) :: f(:), jf(:, :)

    call curv_modine_jacobian_inv(sb_p, cv_p, y, f, jf)
    f(:) = f(:) - x_p(:)
  end subroutine getf 


  ! ---------------------------------------------------------
  subroutine curv_modine_x2chi(sb, cv, x, chi)
    type(simul_box_t),   target, intent(in)  :: sb
    type(curv_modine_t), target, intent(in)  :: cv
    FLOAT,                       intent(in)  :: x(:)    ! x(sb%dim)
    FLOAT,                       intent(out) :: chi(:)  ! chi(sb%dim)

    logical :: conv
    type(root_solver_t) :: rs

    call root_solver_init(rs, sb%dim,  &
      solver_type = ROOT_NEWTON, maxiter = 500, abs_tolerance = CNST(1.0e-10))

    sb_p  => sb
    cv_p  => cv
    SAFE_ALLOCATE(x_p(1:sb%dim))
    x_p(:) = x(:)

    call droot_solver_run(rs, getf, chi, conv, startval = x)

    SAFE_DEALLOCATE_A(x_p)
    nullify(sb_p); nullify(cv_p)

    if(.not.conv) then
      message(1) = "During the construction of the adaptive grid, the Newton-Raphson"
      message(2) = "method did not converge for point:"
      write(message(3),'(3f14.6)') x(1:sb%dim)
      call write_fatal(3)
    end if

  end subroutine curv_modine_x2chi

end module curv_modine_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
