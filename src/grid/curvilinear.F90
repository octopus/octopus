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

module curvilinear_oct_m
  use curv_briggs_oct_m
  use curv_gygi_oct_m
  use curv_modine_oct_m
  use global_oct_m
  use ions_oct_m
  use lalg_adv_oct_m
  use lattice_vectors_oct_m
  use messages_oct_m
  use namespace_oct_m
  use parser_oct_m
  use profiling_oct_m
  use root_solver_oct_m
  use simul_box_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use varinfo_oct_m

  implicit none

  private
  public ::                      &
    curvilinear_t,               &
    curvilinear_init,            &
    curvilinear_copy,            &
    curvilinear_end,             &
    curvilinear_chi2x,           &
    curvilinear_x2chi,           &
    curvilinear_det_Jac,         &
    curvilinear_write_info

  integer, parameter, public :: &
    CURV_METHOD_UNIFORM = 1,    &
    CURV_METHOD_GYGI    = 2,    &
    CURV_METHOD_BRIGGS  = 3,    &
    CURV_METHOD_MODINE  = 4

  type curvilinear_t
    private
    integer, public :: method
    type(curv_gygi_t)   :: gygi
    type(curv_briggs_t) :: briggs
    type(curv_modine_t) :: modine
    type(root_solver_t) :: rs
    FLOAT, public :: min_mesh_scaling_product ! product of the smallest scaling :: min(distance between the grid points / spacing)
  end type curvilinear_t

  character(len=23), parameter :: dump_tag = '*** curvilinear_dump **'

contains

  ! ---------------------------------------------------------
  subroutine curvilinear_init(cv, namespace, sb, ions, spacing)
    type(curvilinear_t), intent(out) :: cv
    type(namespace_t),   intent(in)  :: namespace
    type(simul_box_t),   intent(in)  :: sb
    type(ions_t),        intent(in)  :: ions
    FLOAT,               intent(in)  :: spacing(:)

    PUSH_SUB(curvilinear_init)

    !%Variable CurvMethod
    !%Type integer
    !%Default curv_uniform
    !%Section Mesh::Curvilinear
    !%Description
    !% The relevant functions in octopus are represented on a mesh in real space.
    !% This mesh may be an evenly spaced regular rectangular grid (standard mode),
    !% or else an adaptive or curvilinear grid. We have implemented
    !% three kinds of adaptive meshes, although only one is currently working,
    !% the one invented by F. Gygi (<tt>curv_gygi</tt>). The code will stop if any of
    !% the other two is invoked. All are experimental with domain parallelization.
    !%Option curv_uniform 1
    !% Regular, uniform rectangular grid.
    !%Option curv_gygi 2
    !% The deformation of the grid is done according to the scheme described by
    !% F. Gygi [F. Gygi and G. Galli, <i>Phys. Rev. B</i> <b>52</b>, R2229 (1995)].
    !%Option curv_briggs 3
    !% The deformation of the grid is done according to the scheme described by
    !% Briggs [E.L. Briggs, D.J. Sullivan, and J. Bernholc, <i>Phys. Rev. B</i> <b>54</b> 14362 (1996)]
    !% (NOT WORKING).
    !%Option curv_modine 4
    !% The deformation of the grid is done according to the scheme described by
    !% Modine [N.A. Modine, G. Zumbach and E. Kaxiras, <i>Phys. Rev. B</i> <b>55</b>, 10289 (1997)]
    !% (NOT WORKING).
    !%End
    call parse_variable(namespace, 'CurvMethod', CURV_METHOD_UNIFORM, cv%method)
    if(.not.varinfo_valid_option('CurvMethod', cv%method)) call messages_input_error(namespace, 'CurvMethod')
    call messages_print_var_option(stdout, "CurvMethod", cv%method)

    ! FIXME: The other two methods are apparently not working
    if(cv%method > CURV_METHOD_GYGI) call messages_experimental('Selected curvilinear coordinates method')

    select case(cv%method)
    case(CURV_METHOD_GYGI)
      call curv_gygi_init(cv%gygi, namespace, sb, ions, cv%min_mesh_scaling_product)
    case(CURV_METHOD_BRIGGS)
      call curv_briggs_init(cv%briggs, namespace, sb, spacing, cv%min_mesh_scaling_product)
    case(CURV_METHOD_MODINE)
      call curv_modine_init(cv%modine, namespace, sb, ions, spacing, cv%min_mesh_scaling_product)
    end select

    ! initialize root solver
    call root_solver_init(cv%rs, namespace, sb%dim,  &
      solver_type = ROOT_NEWTON, maxiter = 500, abs_tolerance = CNST(1.0e-10))

    POP_SUB(curvilinear_init)
  end subroutine curvilinear_init

  ! -------------------------------------------------------------- 
  subroutine curvilinear_copy(this_out, this_in)
    type(curvilinear_t), intent(inout) :: this_out
    type(curvilinear_t), intent(in)    :: this_in
    !
    PUSH_SUB(curvilinear_copy)
    this_out%method=this_in%method
    select case(this_in%method)
    case(CURV_METHOD_GYGI)
      call curv_gygi_copy(this_out%gygi, this_in%gygi)
    case(CURV_METHOD_BRIGGS)
      call curv_briggs_copy(this_out%briggs, this_in%briggs)
    case(CURV_METHOD_MODINE)
      call curv_modine_copy(this_out%modine, this_in%modine)
    end select
    POP_SUB(curvilinear_copy)
    return
  end subroutine curvilinear_copy

  ! ---------------------------------------------------------
  subroutine curvilinear_end(cv)
    type(curvilinear_t), intent(inout) :: cv

    PUSH_SUB(curvilinear_end)

    select case(cv%method)
    case(CURV_METHOD_GYGI)
      call curv_gygi_end(cv%gygi)
    case(CURV_METHOD_BRIGGS)
      !
    case(CURV_METHOD_MODINE)
      call curv_modine_end(cv%modine)
    end select

    POP_SUB(curvilinear_end)
  end subroutine curvilinear_end


  ! ---------------------------------------------------------
  subroutine curvilinear_chi2x(sb, latt, cv, chi, x)
    type(simul_box_t),       intent(in)  :: sb
    type(lattice_vectors_t), intent(in)  :: latt
    type(curvilinear_t),     intent(in)  :: cv
    FLOAT,                   intent(in)  :: chi(:)  !< chi(1:sb%dim)
    FLOAT,                   intent(out) :: x(:)    !< x(1:sb%dim)

    ! no push_sub because called too frequently
    x = M_ZERO

    select case(cv%method)
    case(CURV_METHOD_UNIFORM)
      x(1:sb%dim) = matmul(latt%rlattice_primitive(1:sb%dim,1:sb%dim), chi(1:sb%dim))
    case(CURV_METHOD_GYGI)
      call curv_gygi_chi2x(sb, cv%gygi, cv%rs, chi, x)
    case(CURV_METHOD_BRIGGS)
      call curv_briggs_chi2x(sb, cv%briggs, chi, x)
    case(CURV_METHOD_MODINE)
      call curv_modine_chi2x(sb, cv%modine, chi, x)
    end select

  end subroutine curvilinear_chi2x


  ! ---------------------------------------------------------
  subroutine curvilinear_x2chi(sb, latt, cv, x, chi)
    type(simul_box_t),       intent(in)  :: sb
    type(lattice_vectors_t), intent(in)  :: latt
    type(curvilinear_t),     intent(in)  :: cv
    FLOAT,                   intent(in)  :: x(MAX_DIM)
    FLOAT,                   intent(out) :: chi(MAX_DIM)

    PUSH_SUB(curvilinear_x2chi)

    chi = M_ZERO

    select case(cv%method)
    case(CURV_METHOD_UNIFORM)
      chi(1:sb%dim) = matmul(x(1:sb%dim), latt%klattice_primitive)
    case(CURV_METHOD_GYGI)
      call curv_gygi_x2chi(sb, cv%gygi, x, chi)
    case(CURV_METHOD_BRIGGS, CURV_METHOD_MODINE)
      message(1) = "Internal error in curvilinear_x2chi"
      call messages_fatal(1)
    end select

    POP_SUB(curvilinear_x2chi)
  end subroutine curvilinear_x2chi


  ! ---------------------------------------------------------
  FLOAT function curvilinear_det_Jac(sb, latt, cv, x, chi) result(jdet)
    type(simul_box_t),       intent(in)  :: sb
    type(lattice_vectors_t), intent(in)  :: latt
    type(curvilinear_t),     intent(in)  :: cv
    FLOAT,                   intent(in)  :: x(:)    !<   x(sb%dim)
    FLOAT,                   intent(in)  :: chi(:)  !< chi(sb%dim)

    FLOAT :: dummy(MAX_DIM)
    FLOAT, allocatable :: Jac(:,:)
    integer :: i

    ! No PUSH_SUB, called too often

    SAFE_ALLOCATE(Jac(1:sb%dim, 1:sb%dim))

    select case(cv%method)
    case(CURV_METHOD_UNIFORM)
      Jac(1:sb%dim, 1:sb%dim) = latt%rlattice_primitive(1:sb%dim, 1:sb%dim)
      jdet = lalg_determinant(sb%dim, Jac, preserve_mat = .false.)      
    case(CURV_METHOD_GYGI)
      call curv_gygi_jacobian(sb, cv%gygi, x, dummy, Jac)
      jdet = M_ONE/lalg_determinant(sb%dim, Jac, preserve_mat = .false.)
    case(CURV_METHOD_BRIGGS)
      call curv_briggs_jacobian_inv(sb, cv%briggs, chi, Jac)
      jdet = M_ONE
      do i = 1, sb%dim
        jdet = jdet * Jac(i,i) ! Jacobian is diagonal in this method
      end do
    case(CURV_METHOD_MODINE)
      call curv_modine_jacobian_inv(sb, cv%modine, chi, dummy, Jac)
      jdet = M_ONE*lalg_determinant(sb%dim, Jac, preserve_mat = .false.)
    end select

    SAFE_DEALLOCATE_A(Jac)

  end function curvilinear_det_Jac

  ! ---------------------------------------------------------
  subroutine curvilinear_write_info(cv, unit)
    type(curvilinear_t), intent(in) :: cv
    integer,            intent(in) :: unit

    PUSH_SUB(curvilinear_write_info)

    select case(cv%method)
    case(CURV_METHOD_GYGI)
      write(message(1), '(a)')  '  Curvilinear Method = gygi'
      write(message(2), '(a)')  '  Gygi Parameters:'
      write(message(3), '(4x,a,f6.3)')  'A = ', cv%gygi%a
      write(message(4), '(4x,3a,f6.3)') 'alpha [', &
        trim(units_abbrev(units_out%length)), '] = ', &
        units_from_atomic(units_out%length, cv%gygi%alpha)
      write(message(5), '(4x,3a,f6.3)') 'beta  [', &
        trim(units_abbrev(units_out%length)), '] = ', &
        units_from_atomic(units_out%length, cv%gygi%beta)
      call messages_info(5, unit)

    case(CURV_METHOD_BRIGGS)
      write(message(1), '(a)') '  Curvilinear Method = briggs'
      call messages_info(1, unit)

    case(CURV_METHOD_MODINE)
      write(message(1), '(a)') ' Curvilinear  Method = modine'
      call messages_info(1, unit)

    end select

    POP_SUB(curvilinear_write_info)
  end subroutine curvilinear_write_info

  ! ---------------------------------------------------------

end module curvilinear_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
