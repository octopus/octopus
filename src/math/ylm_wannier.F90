  !              
  ! Copyright (C) 2003-2013 Quantum ESPRESSO and Wannier90 groups
  ! This file is distributed under the terms of the
  ! GNU General Public License. See the file `License'
  ! in the root directory of the present distribution,
  ! or http://www.gnu.org/copyleft/gpl.txt .
  !
  ! routine adapted from the original Pwscf interface PP/pw2wannier90.f90
  !
#include "global.h"

module ylm_wannier_oct_m
  use global_oct_m
  use loct_math_oct_m


  implicit none

  private

  public :: ylm_wannier

contains
  subroutine ylm_wannier(ylm,l,mr,r,nr)
  ! this routine returns in ylm(r) the values at the nr points r(1:3,1:nr)
  ! of the spherical harmonic identified  by indices (l,mr)
  ! in table 3.1 of the wannierf90 specification.
  !
  ! No reference to the particular ylm ordering internal to octopus
  ! is assumed.
  !
  ! If ordering in wannier90 code is changed or extended this should be the
  ! only place to be modified accordingly
  !
  
  ! I/O variables
  !
     INTEGER :: l, mr, nr
     FLOAT :: ylm(nr), r(3,nr)
  !
  ! local variables
  !
  !     FLOAT, EXTERNAL ::  s, p_z,px,py, dz2, dxz, dyz, dx2my2, dxy 
  !     FLOAT, EXTERNAL :: fz3, fxz2, fyz2, fzx2my2, fxyz, fxx2m3y2, fy3x2my2
     FLOAT :: rr, cost, phi
     INTEGER :: ir
     FLOAT :: bs2, bs3, bs6, bs12
     bs2 = M_ONE/sqrt(M_TWO)
     bs3=M_ONE/sqrt(CNST(3.0))
     bs6 = M_ONE/sqrt(6.d0)
     bs12 = M_ONE/sqrt(CNST(12.0))
  
     IF (l > 3 .or. l < -5 ) stop 'ylm_wannier l out of range '
     IF (l>=0) THEN
        IF (mr < 1 .or. mr > 2*l+1) stop 'ylm_wannier mr out of range'
     ELSE
        IF (mr < 1 .or. mr > abs(l)+1 ) stop 'ylm_wannier mr out of range'
     endIF
  
     DO ir=1, nr
        rr = sqrt( r(1,ir)*r(1,ir) +  r(2,ir)*r(2,ir) + r(3,ir)*r(3,ir) )
        IF (rr < M_EPSILON) then !stop 'ylm_wannier rr too small'
           ylm(ir)=0
           cycle
        end IF
        cost =  r(3,ir) / rr
        !
        !  beware the arc tan, it is defined modulo M_PI
        !
        IF (r(1,ir) > M_EPSILON) THEN
           phi = atan( r(2,ir)/r(1,ir) )
        ELSEIF (r(1,ir) < -M_EPSILON ) THEN
           phi = atan( r(2,ir)/r(1,ir) ) + M_PI
        ELSE
           phi = sign( M_PI/M_TWO,r(2,ir) )
        endIF
  
  
        IF (l==0) THEN   ! s orbital
           ylm(ir) = s(cost,phi)
        endIF
        IF (l==1) THEN   ! p orbitals
           IF (mr==1) ylm(ir) = p_z(cost,phi)
           IF (mr==2) ylm(ir) = px(cost,phi)
           IF (mr==3) ylm(ir) = py(cost,phi)
        endIF
        IF (l==2) THEN   ! d orbitals
           IF (mr==1) ylm(ir) = dz2(cost,phi)
           IF (mr==2) ylm(ir) = dxz(cost,phi)
           IF (mr==3) ylm(ir) = dyz(cost,phi)
           IF (mr==4) ylm(ir) = dx2my2(cost,phi)
           IF (mr==5) ylm(ir) = dxy(cost,phi)
        endIF
        IF (l==3) THEN   ! f orbitals
           IF (mr==1) ylm(ir) = fz3(cost,phi)
           IF (mr==2) ylm(ir) = fxz2(cost,phi)
           IF (mr==3) ylm(ir) = fyz2(cost,phi)
           IF (mr==4) ylm(ir) = fzx2my2(cost,phi)
           IF (mr==5) ylm(ir) = fxyz(cost,phi)
           IF (mr==6) ylm(ir) = fxx2m3y2(cost,phi)
           IF (mr==7) ylm(ir) = fy3x2my2(cost,phi)
        endIF
        IF (l==-1) THEN  !  sp hybrids
           IF (mr==1) ylm(ir) = bs2 * ( s(cost,phi) + px(cost,phi) )
           IF (mr==2) ylm(ir) = bs2 * ( s(cost,phi) - px(cost,phi) )
        endIF
        IF (l==-2) THEN  !  sp2 hybrids
           IF (mr==1) ylm(ir) = bs3*s(cost,phi)-bs6*px(cost,phi)+bs2*py(cost,phi)
           IF (mr==2) ylm(ir) = bs3*s(cost,phi)-bs6*px(cost,phi)-bs2*py(cost,phi)
           IF (mr==3) ylm(ir) = bs3*s(cost,phi) +M_TWO*bs6*px(cost,phi)
        endIF
        IF (l==-3) THEN  !  sp3 hybrids
           IF (mr==1) ylm(ir) = 0.5d0*(s(cost,phi)+px(cost,phi)+py(cost,phi)+p_z(cost,phi))
           IF (mr==2) ylm(ir) = 0.5d0*(s(cost,phi)+px(cost,phi)-py(cost,phi)-p_z(cost,phi))
           IF (mr==3) ylm(ir) = 0.5d0*(s(cost,phi)-px(cost,phi)+py(cost,phi)-p_z(cost,phi))
           IF (mr==4) ylm(ir) = 0.5d0*(s(cost,phi)-px(cost,phi)-py(cost,phi)+p_z(cost,phi))
        endIF
        IF (l==-4) THEN  !  sp3d hybrids
           IF (mr==1) ylm(ir) = bs3*s(cost,phi)-bs6*px(cost,phi)+bs2*py(cost,phi)
           IF (mr==2) ylm(ir) = bs3*s(cost,phi)-bs6*px(cost,phi)-bs2*py(cost,phi)
           IF (mr==3) ylm(ir) = bs3*s(cost,phi) +M_TWO*bs6*px(cost,phi)
           IF (mr==4) ylm(ir) = bs2*p_z(cost,phi)+bs2*dz2(cost,phi)
           IF (mr==5) ylm(ir) =-bs2*p_z(cost,phi)+bs2*dz2(cost,phi)
        endIF
        IF (l==-5) THEN  ! sp3d2 hybrids
           IF (mr==1) ylm(ir) = bs6*s(cost,phi)-bs2*px(cost,phi)-bs12*dz2(cost,phi)+.5d0*dx2my2(cost,phi)
           IF (mr==2) ylm(ir) = bs6*s(cost,phi)+bs2*px(cost,phi)-bs12*dz2(cost,phi)+.5d0*dx2my2(cost,phi)
           IF (mr==3) ylm(ir) = bs6*s(cost,phi)-bs2*py(cost,phi)-bs12*dz2(cost,phi)-.5d0*dx2my2(cost,phi)
           IF (mr==4) ylm(ir) = bs6*s(cost,phi)+bs2*py(cost,phi)-bs12*dz2(cost,phi)-.5d0*dx2my2(cost,phi)
           IF (mr==5) ylm(ir) = bs6*s(cost,phi)-bs2*p_z(cost,phi)+bs3*dz2(cost,phi)
           IF (mr==6) ylm(ir) = bs6*s(cost,phi)+bs2*p_z(cost,phi)+bs3*dz2(cost,phi)
        endIF
  
     endDO
  
     return
   end subroutine ylm_wannier
  
    !======== l = 0 =====================================================================
     function s(cost,phi)
       FLOAT :: s, cost,phi
       s = M_ONE/ sqrt((M_FOUR*M_PI))
       return
    end function s

    !======== l = 1 =====================================================================
    function p_z(cost,phi)
       FLOAT ::p_z, cost,phi
       p_z =  sqrt(CNST(3.0)/(M_FOUR*M_PI)) * cost
       return
    end function p_z

    function px(cost,phi)
       FLOAT ::px, cost, phi, sint
       sint = sqrt(abs(M_ONE - cost*cost))
       px =  sqrt(CNST(3.0)/(M_FOUR*M_PI)) * sint * cos(phi)
       return
    end function px

    function py(cost,phi)
       FLOAT ::py, cost, phi, sint
       sint = sqrt(abs(M_ONE - cost*cost))
       py =  sqrt(CNST(3.0)/(M_FOUR*M_PI)) * sint * sin(phi)
       return
    end function py

    !======== l = 2 =====================================================================
    function dz2(cost,phi)
       FLOAT ::dz2, cost, phi
       dz2 =  sqrt(1.25d0/(M_FOUR*M_PI)) * (CNST(3.0)* cost*cost-M_ONE)
       return
    end function dz2

    function dxz(cost,phi)
       FLOAT ::dxz, cost, phi, sint
       sint = sqrt(abs(M_ONE - cost*cost))
       dxz =  sqrt(15.d0/(M_FOUR*M_PI)) * sint*cost * cos(phi)
       return
    end function dxz

    function dyz(cost,phi)
       FLOAT ::dyz, cost, phi, sint
       sint = sqrt(abs(M_ONE - cost*cost))
       dyz =  sqrt(15.d0/(M_FOUR*M_PI)) * sint*cost * sin(phi)
       return
    end function dyz

    function dx2my2(cost,phi)
       FLOAT ::dx2my2, cost, phi, sint
       sint = sqrt(abs(M_ONE - cost*cost))
       dx2my2 =  sqrt(3.75d0/(M_FOUR*M_PI)) * sint*sint * cos(M_TWO*phi)
       return
    end function dx2my2

    function dxy(cost,phi)
       FLOAT ::dxy, cost, phi, sint
       sint = sqrt(abs(M_ONE - cost*cost))
       dxy =  sqrt(3.75d0/(M_FOUR*M_PI)) * sint*sint * sin(M_TWO*phi)
       return
    end function dxy

    !======== l = 3 =====================================================================
    function fz3(cost,phi)
       FLOAT ::fz3, cost, phi
       fz3 =  0.25d0*sqrt(7.d0/M_PI) * ( 5.d0 * cost * cost - CNST(3.0) ) * cost
       return
    end function fz3

    function fxz2(cost,phi)
       FLOAT ::fxz2, cost, phi, sint
       sint = sqrt(abs(M_ONE - cost*cost))
       fxz2 =  0.25d0*sqrt(10.5d0/M_PI) * ( 5.d0 * cost * cost - M_ONE ) * sint * cos(phi)
       return
    end function fxz2

    function fyz2(cost,phi)
       FLOAT ::fyz2, cost, phi, sint
       sint = sqrt(abs(M_ONE - cost*cost))
       fyz2 =  0.25d0*sqrt(10.5d0/M_PI) * ( 5.d0 * cost * cost - M_ONE ) * sint * sin(phi)
       return
    end function fyz2

    function fzx2my2(cost,phi)
       FLOAT ::fzx2my2, cost, phi, sint
       sint = sqrt(abs(M_ONE - cost*cost))
       fzx2my2 =  0.25d0*sqrt(105d0/M_PI) * sint * sint * cost * cos(M_TWO*phi)
       return
    end function fzx2my2

    function fxyz(cost,phi)
       FLOAT ::fxyz, cost, phi, sint
       sint = sqrt(abs(M_ONE - cost*cost))
       fxyz =  0.25d0*sqrt(105d0/M_PI) * sint * sint * cost * sin(M_TWO*phi)
       return
    end function fxyz

    function fxx2m3y2(cost,phi)
       FLOAT ::fxx2m3y2, cost, phi, sint
       sint = sqrt(abs(M_ONE - cost*cost))
       fxx2m3y2 =  0.25d0*sqrt(17.5d0/M_PI) * sint * sint * sint * cos(CNST(3.0)*phi)
       return
    end function fxx2m3y2

    function fy3x2my2(cost,phi)
       FLOAT ::fy3x2my2, cost, phi, sint
       sint = sqrt(abs(M_ONE - cost*cost))
       fy3x2my2 =  0.25d0*sqrt(17.5d0/M_PI) * sint * sint * sint * sin(CNST(3.0)*phi)
       return
    end function fy3x2my2

end module ylm_wannier_oct_m
