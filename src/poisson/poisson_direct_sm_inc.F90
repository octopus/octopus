!! Copyright (C) 2002-2011 M. Marques, A. Castro, A. Rubio, G. Bertsch, M. Oliveira
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

!-----------------------------------------------------------------
subroutine dpoisson_solve_direct_sm(this, sm, pot, rho)
  type(poisson_t), intent(in)  :: this
  type(submesh_t), intent(in)  :: sm 
  FLOAT,           intent(out) :: pot(:)
  FLOAT,           intent(in)  :: rho(:)

  FLOAT                :: prefactor
  FLOAT                :: aa1, aa2, aa3, aa4
  integer              :: ip, jp, dim, nthreads
  FLOAT                :: xx1(1:MAX_DIM), xx2(1:MAX_DIM), xx3(1:MAX_DIM), xx4(1:MAX_DIM)

  PUSH_SUB(dpoisson_solve_direct_sm)

  nthreads = 1
#ifdef HAVE_OPENMP
  !$omp parallel
  !$omp master
  nthreads = omp_get_num_threads()
  !$omp end master
  !$omp end parallel
#endif

  ASSERT(.not.this%der%mesh%parallel_in_domains)

  dim = sm%mesh%sb%dim

  select case(dim)
  case(3)
    prefactor = M_TWO*M_PI*(M_THREE/(M_PI*M_FOUR))**(M_TWOTHIRD)
  case(2)
    prefactor = M_TWO*sqrt(M_PI)
  case default
    message(1) = "Internal error: poisson_solve_direct can only be called for 2D or 3D."
    ! why not? all that is needed is the appropriate prefactors to be defined above, actually. then 1D, 4D etc. can be done
    call messages_fatal(1)
  end select

  if(.not. sm%mesh%use_curvilinear) then
    prefactor = prefactor / (sm%mesh%volume_element**(M_ONE/sm%mesh%sb%dim))
  end if

  do ip = 1, sm%np - 4 + 1, 4

    xx1(1:dim) = sm%x(ip    , 1:dim)
    xx2(1:dim) = sm%x(ip + 1, 1:dim)
    xx3(1:dim) = sm%x(ip + 2, 1:dim)
    xx4(1:dim) = sm%x(ip + 3, 1:dim)
      
    if(this%der%mesh%use_curvilinear) then

      aa1 = prefactor*rho(ip    )*sm%mesh%vol_pp(ip    )**(M_ONE - M_ONE/sm%mesh%sb%dim)
      aa2 = prefactor*rho(ip + 1)*sm%mesh%vol_pp(ip + 1)**(M_ONE - M_ONE/sm%mesh%sb%dim)
      aa3 = prefactor*rho(ip + 2)*sm%mesh%vol_pp(ip + 2)**(M_ONE - M_ONE/sm%mesh%sb%dim)
      aa4 = prefactor*rho(ip + 3)*sm%mesh%vol_pp(ip + 3)**(M_ONE - M_ONE/sm%mesh%sb%dim)

      !$omp parallel do reduction(+:aa1,aa2,aa3,aa4) schedule(dynamic,sm%np/nthreads)
      do jp = 1, sm%np
        if(ip     /= jp) aa1 = aa1 + rho(jp)/sqrt(sum((xx1(1:dim) - sm%x(jp, 1:dim))**2))*sm%mesh%vol_pp(sm%map(jp))
        if(ip + 1 /= jp) aa2 = aa2 + rho(jp)/sqrt(sum((xx2(1:dim) - sm%x(jp, 1:dim))**2))*sm%mesh%vol_pp(sm%map(jp))
        if(ip + 2 /= jp) aa3 = aa3 + rho(jp)/sqrt(sum((xx3(1:dim) - sm%x(jp, 1:dim))**2))*sm%mesh%vol_pp(sm%map(jp))
        if(ip + 3 /= jp) aa4 = aa4 + rho(jp)/sqrt(sum((xx4(1:dim) - sm%x(jp, 1:dim))**2))*sm%mesh%vol_pp(sm%map(jp))
      end do

    else

      aa1 = prefactor*rho(ip    )
      aa2 = prefactor*rho(ip + 1)
      aa3 = prefactor*rho(ip + 2)
      aa4 = prefactor*rho(ip + 3)

      !$omp parallel do reduction(+:aa1,aa2,aa3,aa4) schedule(dynamic,sm%np/nthreads)
      do jp = 1, sm%np
        if(ip     /= jp) aa1 = aa1 + rho(jp)/sqrt(sum((xx1(1:dim) - sm%x(jp, 1:dim))**2))
        if(ip + 1 /= jp) aa2 = aa2 + rho(jp)/sqrt(sum((xx2(1:dim) - sm%x(jp, 1:dim))**2))
        if(ip + 2 /= jp) aa3 = aa3 + rho(jp)/sqrt(sum((xx3(1:dim) - sm%x(jp, 1:dim))**2))
        if(ip + 3 /= jp) aa4 = aa4 + rho(jp)/sqrt(sum((xx4(1:dim) - sm%x(jp, 1:dim))**2))
      end do
      
    end if

    pot(ip    ) = sm%mesh%volume_element*aa1
    pot(ip + 1) = sm%mesh%volume_element*aa2
    pot(ip + 2) = sm%mesh%volume_element*aa3
    pot(ip + 3) = sm%mesh%volume_element*aa4
    
  end do
  
  
  do ip = ip, sm%np

    aa1 = M_ZERO

    xx1(1:dim) = sm%x(ip,1:dim)
    if(sm%mesh%use_curvilinear) then
      !$omp parallel do reduction(+:aa1)
      do jp = 1, sm%np
        if(ip == jp) then
          aa1 = aa1 + prefactor*rho(ip)*sm%mesh%vol_pp(sm%map(jp))**(M_ONE - M_ONE/sm%mesh%sb%dim)
        else
          aa1 = aa1 + rho(jp)/sqrt(sum((xx1(1:dim) - sm%x(jp, 1:dim))**2))*sm%mesh%vol_pp(sm%map(jp))
        end if
      end do
    else
      !$omp parallel do reduction(+:aa1)
      do jp = 1, sm%np
        if(ip == jp) then
          aa1 = aa1 + prefactor*rho(ip)
        else
          aa1 = aa1 + rho(jp)/sqrt(sum((xx1(1:dim) - sm%x(jp, 1:dim))**2))
        end if
      end do
    end if

    pot(ip) = sm%mesh%volume_element*aa1
    
  end do
  
  POP_SUB(dpoisson_solve_direct_sm) 
end subroutine dpoisson_solve_direct_sm

  subroutine zpoisson_solve_direct_sm(this, sm, pot, rho)
    type(poisson_t),      intent(in)    :: this
    type(submesh_t),      intent(in)    :: sm
    CMPLX,                intent(inout) :: pot(:)  !< pot(mesh%np)
    CMPLX,                intent(in)    :: rho(:)  !< rho(mesh%np)

    FLOAT, allocatable :: aux1(:), aux2(:)
    type(derivatives_t), pointer :: der

    der => this%der

    PUSH_SUB(zpoisson_solve_direct_sm)

    SAFE_ALLOCATE(aux1(1:sm%np))
    SAFE_ALLOCATE(aux2(1:sm%np))
    ! first the real part
    aux1(1:sm%np) = real(rho(1:sm%np))
    aux2(1:sm%np) = real(pot(1:sm%np))
    call dpoisson_solve_direct_sm(this, sm, aux2, aux1)
    pot(1:sm%np)  = aux2(1:sm%np)

    ! now the imaginary part
    aux1(1:sm%np) = aimag(rho(1:sm%np))
    aux2(1:sm%np) = aimag(pot(1:sm%np))
    call dpoisson_solve_direct_sm(this, sm, aux2, aux1)
    pot(1:sm%np) = pot(1:sm%np) + M_zI*aux2(1:sm%np)

    SAFE_DEALLOCATE_A(aux1)
    SAFE_DEALLOCATE_A(aux2)

    POP_SUB(zpoisson_solve_direct_sm)
  end subroutine zpoisson_solve_direct_sm


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
