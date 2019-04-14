  !              
  ! Copyright (C) 2003-2013 Quantum ESPRESSO and Wannier90 groups
  ! This file is distributed under the terms of the
  ! GNU General Public License. See the file `License'
  ! in the root directory of the present distribution,
  ! or http://www.gnu.org/copyleft/gpl.txt .
  !
  ! routine adapted from the original Pwscf interface PP/pw2wannier90.f90
  !
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
  
     IMPLICIT NONE
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
     bs2 = 1.d0/sqrt(2.d0)
     bs3=1.d0/sqrt(3.d0)
     bs6 = 1.d0/sqrt(6.d0)
     bs12 = 1.d0/sqrt(12.d0)
  
     IF (l > 3 .or. l < -5 ) stop 'ylm_wannier l out of range '
     IF (l>=0) THEN
        IF (mr < 1 .or. mr > 2*l+1) stop 'ylm_wannier mr out of range'
     ELSE
        IF (mr < 1 .or. mr > abs(l)+1 ) stop 'ylm_wannier mr out of range'
     ENDIF
  
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
           phi = sign( M_PI/2.d0,r(2,ir) )
        ENDIF
  
  
        IF (l==0) THEN   ! s orbital
           ylm(ir) = s(cost,phi)
        ENDIF
        IF (l==1) THEN   ! p orbitals
           IF (mr==1) ylm(ir) = p_z(cost,phi)
           IF (mr==2) ylm(ir) = px(cost,phi)
           IF (mr==3) ylm(ir) = py(cost,phi)
        ENDIF
        IF (l==2) THEN   ! d orbitals
           IF (mr==1) ylm(ir) = dz2(cost,phi)
           IF (mr==2) ylm(ir) = dxz(cost,phi)
           IF (mr==3) ylm(ir) = dyz(cost,phi)
           IF (mr==4) ylm(ir) = dx2my2(cost,phi)
           IF (mr==5) ylm(ir) = dxy(cost,phi)
        ENDIF
        IF (l==3) THEN   ! f orbitals
           IF (mr==1) ylm(ir) = fz3(cost,phi)
           IF (mr==2) ylm(ir) = fxz2(cost,phi)
           IF (mr==3) ylm(ir) = fyz2(cost,phi)
           IF (mr==4) ylm(ir) = fzx2my2(cost,phi)
           IF (mr==5) ylm(ir) = fxyz(cost,phi)
           IF (mr==6) ylm(ir) = fxx2m3y2(cost,phi)
           IF (mr==7) ylm(ir) = fy3x2my2(cost,phi)
        ENDIF
        IF (l==-1) THEN  !  sp hybrids
           IF (mr==1) ylm(ir) = bs2 * ( s(cost,phi) + px(cost,phi) )
           IF (mr==2) ylm(ir) = bs2 * ( s(cost,phi) - px(cost,phi) )
        ENDIF
        IF (l==-2) THEN  !  sp2 hybrids
           IF (mr==1) ylm(ir) = bs3*s(cost,phi)-bs6*px(cost,phi)+bs2*py(cost,phi)
           IF (mr==2) ylm(ir) = bs3*s(cost,phi)-bs6*px(cost,phi)-bs2*py(cost,phi)
           IF (mr==3) ylm(ir) = bs3*s(cost,phi) +2.d0*bs6*px(cost,phi)
        ENDIF
        IF (l==-3) THEN  !  sp3 hybrids
           IF (mr==1) ylm(ir) = 0.5d0*(s(cost,phi)+px(cost,phi)+py(cost,phi)+p_z(cost,phi))
           IF (mr==2) ylm(ir) = 0.5d0*(s(cost,phi)+px(cost,phi)-py(cost,phi)-p_z(cost,phi))
           IF (mr==3) ylm(ir) = 0.5d0*(s(cost,phi)-px(cost,phi)+py(cost,phi)-p_z(cost,phi))
           IF (mr==4) ylm(ir) = 0.5d0*(s(cost,phi)-px(cost,phi)-py(cost,phi)+p_z(cost,phi))
        ENDIF
        IF (l==-4) THEN  !  sp3d hybrids
           IF (mr==1) ylm(ir) = bs3*s(cost,phi)-bs6*px(cost,phi)+bs2*py(cost,phi)
           IF (mr==2) ylm(ir) = bs3*s(cost,phi)-bs6*px(cost,phi)-bs2*py(cost,phi)
           IF (mr==3) ylm(ir) = bs3*s(cost,phi) +2.d0*bs6*px(cost,phi)
           IF (mr==4) ylm(ir) = bs2*p_z(cost,phi)+bs2*dz2(cost,phi)
           IF (mr==5) ylm(ir) =-bs2*p_z(cost,phi)+bs2*dz2(cost,phi)
        ENDIF
        IF (l==-5) THEN  ! sp3d2 hybrids
           IF (mr==1) ylm(ir) = bs6*s(cost,phi)-bs2*px(cost,phi)-bs12*dz2(cost,phi)+.5d0*dx2my2(cost,phi)
           IF (mr==2) ylm(ir) = bs6*s(cost,phi)+bs2*px(cost,phi)-bs12*dz2(cost,phi)+.5d0*dx2my2(cost,phi)
           IF (mr==3) ylm(ir) = bs6*s(cost,phi)-bs2*py(cost,phi)-bs12*dz2(cost,phi)-.5d0*dx2my2(cost,phi)
           IF (mr==4) ylm(ir) = bs6*s(cost,phi)+bs2*py(cost,phi)-bs12*dz2(cost,phi)-.5d0*dx2my2(cost,phi)
           IF (mr==5) ylm(ir) = bs6*s(cost,phi)-bs2*p_z(cost,phi)+bs3*dz2(cost,phi)
           IF (mr==6) ylm(ir) = bs6*s(cost,phi)+bs2*p_z(cost,phi)+bs3*dz2(cost,phi)
        ENDIF
  
     ENDDO
  
     RETURN
   end subroutine ylm_wannier
  
    !======== l = 0 =====================================================================
     FUNCTION s(cost,phi)
       IMPLICIT NONE
       FLOAT :: s, cost,phi
       s = 1.d0/ sqrt((M_FOUR*M_PI))
       RETURN
    END FUNCTION s
    !======== l = 1 =====================================================================
    FUNCTION p_z(cost,phi)
       IMPLICIT NONE
       FLOAT ::p_z, cost,phi
       p_z =  sqrt(3.d0/(M_FOUR*M_PI)) * cost
       RETURN
    END FUNCTION p_z
    FUNCTION px(cost,phi)
       IMPLICIT NONE
       FLOAT ::px, cost, phi, sint
       sint = sqrt(abs(1.d0 - cost*cost))
       px =  sqrt(3.d0/(M_FOUR*M_PI)) * sint * cos(phi)
       RETURN
    END FUNCTION px
    FUNCTION py(cost,phi)
       IMPLICIT NONE
       FLOAT ::py, cost, phi, sint
       sint = sqrt(abs(1.d0 - cost*cost))
       py =  sqrt(3.d0/(M_FOUR*M_PI)) * sint * sin(phi)
       RETURN
    END FUNCTION py
    !======== l = 2 =====================================================================
    FUNCTION dz2(cost,phi)
       IMPLICIT NONE
       FLOAT ::dz2, cost, phi
       dz2 =  sqrt(1.25d0/(M_FOUR*M_PI)) * (3.d0* cost*cost-1.d0)
       RETURN
    END FUNCTION dz2
    FUNCTION dxz(cost,phi)
       IMPLICIT NONE
       FLOAT ::dxz, cost, phi, sint
       sint = sqrt(abs(1.d0 - cost*cost))
       dxz =  sqrt(15.d0/(M_FOUR*M_PI)) * sint*cost * cos(phi)
       RETURN
    END FUNCTION dxz
    FUNCTION dyz(cost,phi)
       IMPLICIT NONE
       FLOAT ::dyz, cost, phi, sint
       sint = sqrt(abs(1.d0 - cost*cost))
       dyz =  sqrt(15.d0/(M_FOUR*M_PI)) * sint*cost * sin(phi)
       RETURN
    END FUNCTION dyz
    FUNCTION dx2my2(cost,phi)
       IMPLICIT NONE
       FLOAT ::dx2my2, cost, phi, sint
       sint = sqrt(abs(1.d0 - cost*cost))
       dx2my2 =  sqrt(3.75d0/(M_FOUR*M_PI)) * sint*sint * cos(2.d0*phi)
       RETURN
    END FUNCTION dx2my2
    FUNCTION dxy(cost,phi)
       IMPLICIT NONE
       FLOAT ::dxy, cost, phi, sint
       sint = sqrt(abs(1.d0 - cost*cost))
       dxy =  sqrt(3.75d0/(M_FOUR*M_PI)) * sint*sint * sin(2.d0*phi)
       RETURN
    END FUNCTION dxy
    !======== l = 3 =====================================================================
    FUNCTION fz3(cost,phi)
       IMPLICIT NONE
       FLOAT ::fz3, cost, phi
       fz3 =  0.25d0*sqrt(7.d0/M_PI) * ( 5.d0 * cost * cost - 3.d0 ) * cost
       RETURN
    END FUNCTION fz3
    FUNCTION fxz2(cost,phi)
       IMPLICIT NONE
       FLOAT ::fxz2, cost, phi, sint
       sint = sqrt(abs(1.d0 - cost*cost))
       fxz2 =  0.25d0*sqrt(10.5d0/M_PI) * ( 5.d0 * cost * cost - 1.d0 ) * sint * cos(phi)
       RETURN
    END FUNCTION fxz2
    FUNCTION fyz2(cost,phi)
       IMPLICIT NONE
       FLOAT ::fyz2, cost, phi, sint
       sint = sqrt(abs(1.d0 - cost*cost))
       fyz2 =  0.25d0*sqrt(10.5d0/M_PI) * ( 5.d0 * cost * cost - 1.d0 ) * sint * sin(phi)
       RETURN
    END FUNCTION fyz2
    FUNCTION fzx2my2(cost,phi)
       IMPLICIT NONE
       FLOAT ::fzx2my2, cost, phi, sint
       sint = sqrt(abs(1.d0 - cost*cost))
       fzx2my2 =  0.25d0*sqrt(105d0/M_PI) * sint * sint * cost * cos(2.d0*phi)
       RETURN
    END FUNCTION fzx2my2
    FUNCTION fxyz(cost,phi)
       IMPLICIT NONE
       FLOAT ::fxyz, cost, phi, sint
       sint = sqrt(abs(1.d0 - cost*cost))
       fxyz =  0.25d0*sqrt(105d0/M_PI) * sint * sint * cost * sin(2.d0*phi)
       RETURN
    END FUNCTION fxyz
    FUNCTION fxx2m3y2(cost,phi)
       IMPLICIT NONE
       FLOAT ::fxx2m3y2, cost, phi, sint
       sint = sqrt(abs(1.d0 - cost*cost))
       fxx2m3y2 =  0.25d0*sqrt(17.5d0/M_PI) * sint * sint * sint * cos(3.d0*phi)
       RETURN
    END FUNCTION fxx2m3y2
    FUNCTION fy3x2my2(cost,phi)
       IMPLICIT NONE
       FLOAT ::fy3x2my2, cost, phi, sint
       sint = sqrt(abs(1.d0 - cost*cost))
       fy3x2my2 =  0.25d0*sqrt(17.5d0/M_PI) * sint * sint * sint * sin(3.d0*phi)
       RETURN
    END FUNCTION fy3x2my2

