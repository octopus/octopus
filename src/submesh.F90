!! Copyright (C) 2007 X. Andrade
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
!! $Id: submesh.F90 2781 2007-03-23 10:58:32Z lorenzen $
 
#include "global.h"
  
module submesh_m
  use global_m
  use messages_m
  use mesh_m
  use mpi_m
  use simul_box_m

  implicit none
  private 

  public ::                 &
       submesh_t,           &
       submesh_null,        &
       submesh_init_sphere, &
       submesh_copy,        &
       dsm_integrate,       &
       zsm_integrate,       &
       ddsm_integrate_prod,  &
       zdsm_integrate_prod,  &
       zzsm_integrate_prod,  &
       submesh_end
 
  type submesh_t
     integer          :: ns = 0        ! number of points inside the submesh
     integer          :: np_part
     integer, pointer :: jxyz(:)        ! index in the mesh of the points inside the sphere
     integer, pointer :: jxyz_inv(:)    ! and the inverse
  end type submesh_t
  
contains
  
  subroutine submesh_null(sm)
    type(submesh_t), intent(out) :: sm

    sm%ns = 0
    sm%np_part = 0
    nullify(sm%jxyz)
    nullify(sm%jxyz_inv)

  end subroutine submesh_null

  subroutine submesh_init_sphere(this, sb, m, center, rc)
    type(submesh_t),   intent(out)  :: this
    type(simul_box_t), intent(in)   :: sb
    type(mesh_t),      intent(in)   :: m
    FLOAT,             intent(in)   :: center(1:MAX_DIM)
    FLOAT,             intent(in)   :: rc
    
    FLOAT :: r
    integer :: ii, is, kk, ip

    ASSERT(this%ns == 0)

    call push_sub('submesh.submesh_init_sphere')

    this%np_part = m%np_part

    if (any(rc + m%h(1) >= sb%lsize(1:sb%periodic_dim))) then
      message(1)='Sphere is larger than the box size'
      write(message(2),'(a,f12.6,a)')  '  rc_max+h = ', rc + m%h(1), ' [b]'
      write(message(3),'(a,3f12.4,a)') '  lsize    = ', sb%lsize, ' [b]'
      message(4)='Please change pseudopotential'
      call write_fatal(4)
    end if

    ! Get the total number of points inside the sphere
    is = 0
    do ip = 1, m%np
      do ii = 1, 3**sb%periodic_dim
        call mesh_r(m, ip, r, a=center(:) + sb%shift(ii,:))
        if(r > rc + m%h(1)) cycle
        is = is + 1
        exit
      end do
    end do

    this%ns = is
    ALLOCATE(this%jxyz(this%ns), this%ns)

    ! Get the index of the points inside the sphere
    is = 0
    do kk = 1, m%np
      do ii = 1, 3**sb%periodic_dim
        call mesh_r(m, kk, r, a=center + sb%shift(ii,:))
        ! we enlarge slightly the mesh (good for the interpolation scheme)
        if(r > rc + m%h(1)) cycle
        is = is + 1
        this%jxyz(is) = kk
        exit
      end do
    end do

    !build the inverse of jxyz, points not in the sphere go to 0
    ALLOCATE(this%jxyz_inv(0:m%np_part), m%np_part+1)

    this%jxyz_inv(0:m%np_part) = 0

    do is = 1, this%ns
      this%jxyz_inv(this%jxyz(is)) = is
    end do
    
    call pop_sub()

  end subroutine submesh_init_sphere

  subroutine submesh_end(this)
    type(submesh_t),   intent(inout)  :: this
    
    call push_sub('submesh.submesh_end')

    if( this%ns /= 0 ) then
      this%ns = 0
      deallocate(this%jxyz)
      deallocate(this%jxyz_inv)
    end if

    call pop_sub()

  end subroutine submesh_end

  subroutine submesh_copy(sm_in, sm_out)
    type(submesh_t),   intent(in)   :: sm_in
    type(submesh_t),   intent(out)  :: sm_out

    call push_sub('submesh.submesh_copy')
    
    ASSERT(sm_out%ns == 0)

    sm_out%ns = sm_in%ns
    sm_out%np_part  = sm_in%np_part
    
    ALLOCATE(sm_out%jxyz(1:sm_out%ns), sm_out%ns)
    ALLOCATE(sm_out%jxyz_inv(0:sm_out%np_part), sm_out%np_part+1)
    
    sm_out%jxyz(1:sm_out%ns)    = sm_in%jxyz(1:sm_in%ns)
    sm_out%jxyz_inv(0:sm_out%np_part) = sm_in%jxyz_inv(0:sm_out%np_part)
   
    call pop_sub()

  end subroutine submesh_copy

  CMPLX function zzsm_integrate_prod(m, sm, f, g) result(res)
    type(mesh_t),    intent(in) :: m
    type(submesh_t), intent(in) :: sm
    CMPLX,           intent(in) :: f(:)
    CMPLX,           intent(in) :: g(:)

#if defined(HAVE_MPI)
    MPI_CMPLX :: tmp
#endif

    res = sum( f(1:sm%ns) * g(1:sm%ns) * m%vol_pp(sm%jxyz(1:sm%ns)) )

#if defined(HAVE_MPI)
    if(m%parallel_in_domains) then
      call MPI_Allreduce(res, tmp, 1, R_MPITYPE, MPI_SUM, m%vp%comm, mpi_err)
      res = tmp
    end if
#endif

  end function zzsm_integrate_prod

#include "undef.F90"
#include "real.F90"
#include "submesh_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "submesh_inc.F90"

end module submesh_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
