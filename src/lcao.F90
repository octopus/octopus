#include "config.h"

module lcao
  use global
  use spline
  use system
  use mesh

  implicit none

contains

!builds a density which is the sum of the atomic densities
subroutine lcao_dens(sys, nspin, rho)
  type(system_type), intent(IN) :: sys
  integer, intent(in) :: nspin
  real(r8), intent(out) :: rho(sys%m%np, nspin)
  
  integer :: ia
  real(r8) :: r
  type(specie_type), pointer :: s
  type(atom_type),   pointer :: a
  
  rho = 0._r8
  do ia = 1, sys%natoms
    a => sys%atom(ia) ! shortcuts
    s => a%spec
    
    select case(s%label(1:5))
    case('jelli', 'point')
      call from_jellium(sys%m, rho(:, 1))
    case default
      call from_pseudopotential(sys%m, rho(:, 1))
    end select
  end do

  ! we now renormalize the density (necessary if we have a charged system)
  ! if spin polarized, we start with paramagnetic density
  r = sys%st%qtot/(dmesh_integrate(sys%m, rho(:, 1))*real(nspin, r8))
  do ia = nspin, 1, -1
    rho(:, ia) = r*rho(:, 1)
  end do

contains
  subroutine from_jellium(m, rho)
    type(mesh_type), intent(in) :: m
    real(r8), intent(inout) :: rho(m%np)

    integer :: i, in_points
    real(r8) :: r

    do i = 1, m%np
      call mesh_r(m, i, r, a=a%x)
      if(r <= s%jradius) then
        in_points = in_points + 1
      end if
    end do
    
    if(in_points > 0) then
      do i = 1, m%np
        call mesh_r(m, i, r, a=a%x)
        if(r <= s%jradius) then
          rho(i) = rho(i) + real(s%Z_val, r8)/(in_points*m%vol_pp)
        end if
      end do
    end if
  end subroutine from_jellium
  
  subroutine from_pseudopotential(m, rho)
    type(mesh_type), intent(in) :: m
    real(r8), intent(inout) :: rho(m%np)

    integer :: i, l
    real(r8) :: r, zel, zval
    R_TYPE :: psi

    do i = 1, m%np
      call mesh_r(m, i, r, a=a%x)
      zval = s%Z_val
      do l = 0 , s%ps%L_max
        if(r >= r_small) then
          psi = splint(s%ps%Ur(l), r)
          zel = min(zval, 2.0_r8*(2*l+1))
          zval = zval - zel
          rho(i) = rho(i) + zel*psi*psi*(r**(2*l))/(4*M_PI)
        end if
      end do
    end do
    
  end subroutine from_pseudopotential
end subroutine lcao_dens

end module lcao
