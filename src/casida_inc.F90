!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
!! $Id: casida_inc.F90 1974 2006-04-03 13:05:44Z acastro $


! ---------------------------------------------------------
function X(ks_matrix_elements) (cas, st, m, dv) result(x)
  type(casida_t), intent(in) :: cas
  type(states_t), intent(in) :: st
  type(mesh_t),   intent(in) :: m
  FLOAT, intent(in)   :: dv(:)
  FLOAT :: x(cas%n_pairs)

  R_TYPE, allocatable :: f(:)
  integer :: k, ia, i, a, sigma

  ALLOCATE(f(m%np), m%np)
  do ia = 1, cas%n_pairs
    i     = cas%pair(ia)%i
    a     = cas%pair(ia)%a
    sigma = cas%pair(ia)%sigma
    do k = 1, m%np
      f(k) = dv(k) * R_CONJ(st%X(psi) (k, 1, i, sigma)) * st%X(psi) (k, 1, a, sigma)
    end do
    x(ia) = X(mf_integrate)(m, f)
  end do

  deallocate(f)
end function X(ks_matrix_elements)

! ---------------------------------------------------------
R_TYPE function X(transition_matrix_element) (cas, ia, x) result(z)
  type(casida_t), intent(in) :: cas
  integer,           intent(in) :: ia
  R_TYPE,            intent(in) :: x(:)

  integer :: jb

  z = R_TOTYPE(M_ZERO)
  if(cas%w(ia) > M_ZERO) then
    do jb = 1, cas%n_pairs
      z = z + x(jb) * (M_ONE/sqrt(cas%s(jb))) * cas%mat(jb, ia)
    end do
    z = (M_ONE/sqrt(cas%w(ia))) * z
  end if

end function X(transition_matrix_element)

! ---------------------------------------------------------
subroutine X(transition_density) (cas, st, m, ia, n0I)
  type(casida_t), intent(in)  :: cas
  type(states_t), intent(in)  :: st
  type(mesh_t),   intent(in)  :: m
  integer,        intent(in)  :: ia
  R_TYPE,         intent(out) :: n0I(:)

  integer :: i, jb
  R_TYPE, allocatable :: x(:)

  call push_sub('casida.transition_density')

  ALLOCATE(x(cas%n_pairs), cas%n_pairs)

  do i = 1, m%np
    do jb = 1, cas%n_pairs
      x(jb) = R_CONJ(st%X(psi)(i, 1, cas%pair(jb)%i, cas%pair(jb)%sigma)) * &
           st%X(psi)(i, 1, cas%pair(jb)%a, cas%pair(jb)%sigma)
    end do
    n0I(i) = X(transition_matrix_element) (cas, ia, x)
  end do

  deallocate(x)
  call pop_sub()
end subroutine X(transition_density)

! ---------------------------------------------------------
subroutine X(get_transition_densities) (cas, sys, trandens)
  type(casida_t),    intent(in) :: cas
  type(system_t),    intent(in) :: sys
  character(len=80), intent(in) :: trandens

  integer :: ia, ierr
  character(len=5) :: intstr
  character(len=130) :: filename
  R_TYPE, allocatable :: n0I(:)
  call push_sub('casida.get_transition_densities')

  ALLOCATE(n0I(sys%gr%m%np), sys%gr%m%np)
  n0I = M_ZERO

  do ia = 1, cas%n_pairs
    if(loct_isinstringlist(ia, trandens)) then
      call X(transition_density) (cas, sys%st, sys%gr%m, ia, n0I)
      write(intstr,'(i5)') ia
      write(intstr,'(i1)') len(trim(adjustl(intstr)))
      write(filename,'(a,i'//trim(intstr)//')') 'n0',ia
      call X(output_function)(sys%outp%how, "linear", trim(filename), &
                              sys%gr%m, sys%gr%sb, n0I, M_ONE, ierr)
    end if
  end do

  deallocate(n0I)
  call pop_sub()
end subroutine X(get_transition_densities)
