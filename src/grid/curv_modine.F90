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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

!> This module implements the curvilinear coordinates given in
!! N. A. Modine, G. Zumbach, and E. Kaxiras, Phys. Rev. B 55, 10289-10301 (1997) 
!!
!! The local refinement was changed for a simple exponential.
!! I believe that the recipe given by the authors is too complicated
!! for me to sort out.

module curv_modine_oct_m
  use global_oct_m
  use ions_oct_m
  use messages_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use root_solver_oct_m
  use simul_box_oct_m
  use unit_oct_m
  use unit_system_oct_m

  implicit none

  private
  public ::                     &
    curv_modine_t,              &
    curv_modine_init,           &
    curv_modine_copy,           &
    curv_modine_end,            &
    curv_modine_chi2x,          &
    curv_modine_jacobian_inv

  type curv_modine_t
    private
    FLOAT              :: L(MAX_DIM) !< size of the box
    FLOAT              :: xbar       !< size of central flat region (in units of L)
    FLOAT              :: Jbar       !< increase in density of points is 1/J

    FLOAT, allocatable :: Jlocal(:)  !< local (around the atoms) refinement
    FLOAT, allocatable :: Jrange(:)  !< local refinement range

    FLOAT, allocatable :: chi_atoms(:,:)
    FLOAT, allocatable :: csi(:,:)

    integer :: natoms
  end type curv_modine_t

  integer, parameter :: qq = 3

  type(simul_box_t),   pointer :: sb_p
  type(curv_modine_t), pointer :: cv_p
  FLOAT,           allocatable :: x_p(:)

contains

  ! ---------------------------------------------------------
  subroutine getf2(csi, ff, jf)
    FLOAT, intent(in)  :: csi(:)
    FLOAT, intent(out) :: ff(:), jf(:, :)

    integer :: i1, j1, i2, j2, index1, index2
    FLOAT :: xx(MAX_DIM), chi2(MAX_DIM), rr, dd, dd2

    PUSH_SUB(getf2)

    ! first we fill in cv%csi with the values we have
    index1 = 1
    do i1 = 1, cv_p%natoms
      do j1 = 1, sb_p%dim
        cv_p%csi(j1, i1) = csi(index1)
        index1 = index1 + 1
      end do
    end do

    ! get ff and jf
    jf(:,:) = M_ZERO
    do i1 = 1, cv_p%natoms
      call curv_modine_chi2chi2(sb_p, cv_p, cv_p%chi_atoms(:,i1), chi2)
      xx(:) = chi2(:)

      do i2 = 1, cv_p%natoms
        rr = sqrt(sum((chi2(:) - cv_p%csi(:,i2))**2))
        dd = exp(-rr**2/(M_TWO*cv_p%Jrange(i2)**2))

        xx(:) = xx(:) - cv_p%Jlocal(i2)*(chi2(:) - cv_p%csi(:,i2)) * dd
      end do

      do j1 = 1, sb_p%dim
        index1 = (i1-1)*sb_p%dim + j1
        ff(index1) = xx(j1) - x_p(index1)

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

    POP_SUB(getf2)
  end subroutine getf2

  ! ---------------------------------------------------------
  subroutine curv_modine_init(cv, namespace, sb, ions, spacing, min_scaling_product)
    type(curv_modine_t), target, intent(out) :: cv
    type(namespace_t),           intent(in)  :: namespace
    type(simul_box_t),   target, intent(in)  :: sb
    type(ions_t),                intent(in)  :: ions
    FLOAT,                       intent(in)  :: spacing(:)
    FLOAT,                       intent(out) :: min_scaling_product

    type(root_solver_t) :: rs

    PUSH_SUB(curv_modine_init)

    !%Variable CurvModineXBar
    !%Type float
    !%Default 1/3
    !%Section Mesh::Curvilinear::Modine
    !%Description
    !% Size of central flat region (in units of <tt>Lsize</tt>). Must be between 0 and 1.
    !% See N. A. Modine, G. Zumbach, and E. Kaxiras, <i>Phys. Rev. B</i> <b>55</b>, 10289-10301 (1997).
    !%End
    call parse_variable(namespace, 'CurvModineXBar', M_ONE/M_THREE, cv%xbar)

    !%Variable CurvModineJBar
    !%Type float
    !%Default 1/2
    !%Section Mesh::Curvilinear::Modine
    !%Description
    !% Increase in density of points is inverse of this parameter.
    !% See N. A. Modine, G. Zumbach, and E. Kaxiras, <i>Phys. Rev. B</i> <b>55</b>, 10289-10301 (1997).
    !%End
    call parse_variable(namespace, 'CurvModineJBar', M_HALF, cv%Jbar)

    cv%L = M_ZERO
    cv%L(1:sb%dim) = sb%lsize(1:sb%dim) / cv%Jbar

    if(cv%xbar<M_ZERO.or.cv%xbar>M_ONE) then
      message(1) = 'The parameter "CurvModineXBar" must lie between 0 and 1.'
      call messages_fatal(1)
    end if

    SAFE_ALLOCATE(cv%Jlocal(1:ions%natoms))
    SAFE_ALLOCATE(cv%Jrange(1:ions%natoms))

    ! \warning: the reading has to be done for each atom kind

    !%Variable CurvModineJlocal
    !%Type float
    !%Default 0.25
    !%Section Mesh::Curvilinear::Modine
    !%Description
    !% Local refinement around the atoms. Must be between 0 and 1.
    !% See N. A. Modine, G. Zumbach, and E. Kaxiras, <i>Phys. Rev. B</i> <b>55</b>, 10289-10301 (1997).
    !%End
    call parse_variable(namespace, 'CurvModineJlocal', CNST(0.25), cv%Jlocal(1))

    !%Variable CurvModineJrange
    !%Type float
    !%Default 2 b
    !%Section Mesh::Curvilinear::Modine
    !%Description
    !% Local refinement range (a length).
    !% See N. A. Modine, G. Zumbach, and E. Kaxiras, <i>Phys. Rev. B</i> <b>55</b>, 10289-10301 (1997).
    !%End
    call parse_variable(namespace, 'CurvModineJrange', M_TWO, cv%Jrange(1), units_inp%length)

    if(cv%Jlocal(1)<M_ZERO.or.cv%Jlocal(1)>M_ONE) then
      message(1) = 'The parameter "CurvModineJlocal" must lie between 0 and 1.'
      call messages_fatal(1)
    end if

    cv%Jlocal(:) = cv%Jlocal(1)
    cv%Jrange(:) = cv%Jrange(1)

    ! initialize root solver for the optimization
    call root_solver_init(rs, namespace, sb%dim, &
        solver_type = ROOT_NEWTON, maxiter = 500, abs_tolerance = CNST(1.0e-10))

    call find_atom_points()
    call optimize()

    cv%natoms = ions%natoms

    call curv_modine_min_scaling(sb, cv, min_scaling_product)

    POP_SUB(curv_modine_init)

  contains

    subroutine find_atom_points()
      integer :: iatom, jj

      PUSH_SUB(curv_modine_init.find_atom_points)

      ! Initialize csi
      SAFE_ALLOCATE(cv%csi(1:sb%dim, 1:ions%natoms))
      cv%csi = ions%pos

      ! get first estimate for chi_atoms
      SAFE_ALLOCATE(cv%chi_atoms(1:sb%dim, 1:ions%natoms))
      do jj = 1, 10  ! \warning: make something better
        do iatom = 1, ions%natoms
          call curv_modine_x2chi(sb, cv, rs, ions%pos(:, iatom), cv%chi_atoms(:, iatom))
        end do
        cv%csi(:,:) = cv%chi_atoms(:,:)
      end do

      do iatom = 1, ions%natoms
        ! These are the chi positions where we want the atoms.
        cv%chi_atoms(:, iatom) = nint(cv%chi_atoms(:, iatom) / spacing(:)) * spacing(:)
      end do

      POP_SUB(curv_modine_init.find_atom_points)
    end subroutine find_atom_points

    subroutine optimize()
      logical :: conv
      integer :: iatom, idim, index
      FLOAT, allocatable :: my_csi(:), start_csi(:)

      PUSH_SUB(curv_modine_init.optimize)

      sb_p  => sb
      cv_p  => cv

      SAFE_ALLOCATE(x_p(1:sb%dim*ions%natoms))
      SAFE_ALLOCATE(my_csi(1:sb%dim*ions%natoms))
      SAFE_ALLOCATE(start_csi(1:sb%dim*ions%natoms))

      do iatom = 1, ions%natoms
        do idim = 1, sb%dim
          index = (iatom-1)*sb%dim + idim
          x_p(index)       = ions%pos(idim, iatom)
          start_csi(index) = cv%chi_atoms(idim, iatom)
        end do
      end do

      call droot_solver_run(rs, getf2, my_csi, conv, startval=start_csi)

      if(.not.conv) then
        message(1) = "During the construction of the adaptive grid, the Newton-Raphson"
        message(2) = "method did not converge."
        call messages_fatal(2)
      end if

      ! Now set csi to the new values
      do iatom = 1, ions%natoms
        do idim = 1, sb%dim
          index = (iatom-1)*sb_p%dim + idim
          cv_p%csi(idim, iatom) = my_csi(index)
        end do
      end do

      SAFE_DEALLOCATE_A(x_p)
      SAFE_DEALLOCATE_A(my_csi)
      SAFE_DEALLOCATE_A(start_csi)

      nullify(sb_p)
      nullify(cv_p)

      POP_SUB(curv_modine_init.optimize)
    end subroutine optimize

  end subroutine curv_modine_init

  ! ---------------------------------------------------------
  subroutine curv_modine_copy(this_out, this_in)
    type(curv_modine_t), intent(inout) :: this_out
    type(curv_modine_t), intent(in)    :: this_in

    PUSH_SUB(curv_modine_copy)

    this_out%L = this_in%L
    this_out%xbar = this_in%xbar
    this_out%Jbar = this_in%Jbar
    SAFE_ALLOCATE_SOURCE_A(this_out%Jlocal, this_in%Jlocal)
    SAFE_ALLOCATE_SOURCE_A(this_out%Jrange, this_in%Jrange)
    SAFE_ALLOCATE_SOURCE_A(this_out%chi_atoms, this_in%chi_atoms)
    SAFE_ALLOCATE_SOURCE_A(this_out%csi, this_in%csi)
    this_out%natoms = this_in%natoms

    POP_SUB(curv_modine_copy)
  end subroutine curv_modine_copy

  ! ---------------------------------------------------------
  subroutine curv_modine_end(cv)
    type(curv_modine_t), intent(inout) :: cv

    PUSH_SUB(curv_modine_end)

    SAFE_DEALLOCATE_A(cv%Jlocal)
    SAFE_DEALLOCATE_A(cv%Jrange)
    SAFE_DEALLOCATE_A(cv%chi_atoms)
    SAFE_DEALLOCATE_A(cv%csi)

    POP_SUB(curv_modine_end)
  end subroutine curv_modine_end


  ! ---------------------------------------------------------
  subroutine curv_modine_chi2chi2(sb, cv, chi_, chi2, Jac)
    type(simul_box_t),   intent(in)  :: sb
    type(curv_modine_t), intent(in)  :: cv
    FLOAT,               intent(in)  :: chi_(:)  !< chi_(sb%dim)
    FLOAT,               intent(out) :: chi2(:)  !< chi2(sb%dim)
    FLOAT,     optional, intent(out) :: Jac(:)   !< the Jacobian of this transformation is diagonal

    FLOAT :: chibar(MAX_DIM), rr, chi
    logical :: neg
    integer :: i

    PUSH_SUB(curv_modine_chi2chi2)

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

    POP_SUB(curv_modine_chi2chi2)
  end subroutine curv_modine_chi2chi2

  ! ---------------------------------------------------------
  subroutine curv_modine_chi2x(sb, cv, chi_, xx)
    type(simul_box_t),   intent(in)  :: sb
    type(curv_modine_t), intent(in)  :: cv
    FLOAT,               intent(in)  :: chi_(:)  !< chi_(sb%dim)
    FLOAT,               intent(out) :: xx(:)    !< xx  (sb%dim)

    FLOAT :: chi2(MAX_DIM), rr, dd
    integer :: iatom

    PUSH_SUB(curv_modine_chi2x)

    call curv_modine_chi2chi2(sb, cv, chi_, chi2)

    xx(:) = chi2(:)
    do iatom = 1, cv%natoms
      rr = max(sqrt(sum((chi2(:) - cv%csi(:, iatom))**2)), CNST(1e-6))
      dd = exp(-rr**2/(M_TWO*cv%Jrange(iatom)**2))

      xx(:) = xx(:) - cv%Jlocal(iatom)*(chi2(:) - cv%csi(:, iatom)) * dd
    end do

    POP_SUB(curv_modine_chi2x)
  end subroutine curv_modine_chi2x


  ! ---------------------------------------------------------
  subroutine curv_modine_jacobian_inv(sb, cv, chi_, xx, Jac)
    type(simul_box_t),   intent(in)  :: sb
    type(curv_modine_t), intent(in)  :: cv
    FLOAT,               intent(in)  :: chi_(:)  !< chi(sb%dim)
    FLOAT,               intent(out) :: xx(:)    !< xx(sb%dim)
    FLOAT,               intent(out) :: Jac(:,:) !< Jac(sb%dim,sb%dim), the Jacobian

    FLOAT :: chi2(MAX_DIM), rr, dd, J2(MAX_DIM)
    integer :: iatom, idim, idim2

    PUSH_SUB(curv_modine_jacobian_inv)

    call curv_modine_chi2chi2(sb, cv, chi_, chi2, J2)

    ! initialize both xx and the Jacobian
    xx(:) = chi2(:)
    Jac(:,:) = M_ZERO
    do idim = 1, sb%dim
      Jac(idim, idim) = M_ONE
    end do

    do iatom = 1, cv%natoms
      rr = max(sqrt(sum(chi2(:) - cv%csi(:, iatom))**2), CNST(1e-6))
      dd = exp(-rr**2/(M_TWO*cv%Jrange(iatom)**2))

      xx(:) = xx(:) -  cv%Jlocal(iatom)*(chi2(:) - cv%csi(:, iatom)) * dd

      do idim = 1, sb%dim
        Jac(idim, idim) = Jac(idim, idim) - cv%Jlocal(iatom) * dd
        do idim2 = 1, sb%dim
          Jac(idim, idim2) = Jac(idim, idim2) + &
             cv%Jlocal(iatom)*(chi2(idim) - cv%csi(idim, iatom))*(chi2(idim2)-cv%csi(idim2, iatom)) * &
             M_TWO/(M_TWO*cv%Jrange(iatom)**2) * dd
        end do
      end do
    end do

    do idim = 1, sb%dim
      Jac(idim, :) = Jac(idim, :) * J2(:)
    end do

    POP_SUB(curv_modine_jacobian_inv)
  end subroutine curv_modine_jacobian_inv


  ! ---------------------------------------------------------
  subroutine getf(yy, ff, jf)
    FLOAT, intent(in)  :: yy(:)
    FLOAT, intent(out) :: ff(:), jf(:, :)

    PUSH_SUB(getf)

    call curv_modine_jacobian_inv(sb_p, cv_p, yy, ff, jf)
    ff(:) = ff(:) - x_p(:)

    POP_SUB(getf)
  end subroutine getf 


  ! ---------------------------------------------------------
  subroutine curv_modine_x2chi(sb, cv, rs, xx, chi)
    type(simul_box_t),   target, intent(in)     :: sb
    type(curv_modine_t), target, intent(in)     :: cv
    type(root_solver_t), target, intent(inout)  :: rs
    FLOAT,                       intent(in)  :: xx(:)   !< xx(sb%dim)
    FLOAT,                       intent(out) :: chi(:)  !< chi(sb%dim)

    logical :: conv

    PUSH_SUB(curv_modine_x2chi)

    sb_p  => sb
    cv_p  => cv
    SAFE_ALLOCATE(x_p(1:sb%dim))
    x_p(:) = xx(:)

    call droot_solver_run(rs, getf, chi, conv, startval = xx)

    SAFE_DEALLOCATE_A(x_p)
    nullify(sb_p)
    nullify(cv_p)

    if(.not.conv) then
      message(1) = "During the construction of the adaptive grid, the Newton-Raphson"
      message(2) = "method did not converge for point:"
      write(message(3),'(3f14.6)') xx(1:sb%dim)
      call messages_fatal(3)
    end if

    POP_SUB(curv_modine_x2chi)
  end subroutine curv_modine_x2chi


  ! ---------------------------------------------------------
  subroutine curv_modine_min_scaling(sb, cv, min_scaling_product)
    type(simul_box_t),   intent(in)  :: sb
    type(curv_modine_t), intent(in)  :: cv
    FLOAT,               intent(out) :: min_scaling_product
    
    min_scaling_product = cv%Jbar**sb%dim
  end subroutine curv_modine_min_scaling

end module curv_modine_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
