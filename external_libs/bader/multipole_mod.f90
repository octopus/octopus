! Copyright 2009
! Wenjie Tang, Andri Arnaldsson, Samuel T. Chill, and Graeme Henkelman
!
! Bader is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! A copy of the GNU General Public License is available at
! http://www.gnu.org/licenses/
!
! The multipole code was writen by:
!
! Sebastien Lebegue <Sebastien.Lebegue@crm2.uhp-nancy.fr>
! Angyan Janos <Janos.Angyan@crm2.uhp-nancy.fr>
! Emmanuel Aubert <emmanuel.aubert@crm2.uhp-nancy.fr>

MODULE multipole_mod

  USE kind_mod
  USE matrix_mod
  USE options_mod
  USE ions_mod
  USE charge_mod
  USE io_mod
  USE chgcar_mod
  USE bader_mod
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: multipole_calc

  CONTAINS

  SUBROUTINE multipole_calc(bdr,ions,chgval,opts)

    TYPE(bader_obj) :: bdr
    TYPE(ions_obj) :: ions
    TYPE(charge_obj) :: chgval
    TYPE(options_obj) :: opts

    REAL(q2),ALLOCATABLE,DIMENSION(:,:) :: dipol
    REAL(q2),ALLOCATABLE,DIMENSION(:) :: charge_check
    REAL(q2),ALLOCATABLE,DIMENSION(:) :: rsquare     
    REAL(q2),ALLOCATABLE,DIMENSION(:,:,:) :: quad
    REAL(q2),DIMENSION(3,3) :: quadtot
    REAL(q2),DIMENSION(3) :: diptot
    REAL(q2) :: chgtot, rsquaretot
    REAL(q2),DIMENSION(3) :: v
    REAL(q2),DIMENSION(3) :: tab
    REAL(q2),DIMENSION(3) :: tab_car
    INTEGER :: n1,n2,n3,i
    INTEGER :: ion_of_point

    REAL(q2),DIMENSION(3) :: grad,voxlen
    REAL(q2) :: rho,vol
    TYPE(charge_obj) :: chgtemp
    TYPE(ions_obj) :: ionstemp
    REAL(q2) :: distance
    REAL(q2),PARAMETER :: charge_sign  = -1.0_q2
    REAL(q2),DIMENSION(3) :: rr0
    REAL(q2),PARAMETER :: bohr2angs  = 0.529177_q2

    IF (opts%ref_flag) THEN
      CALL read_charge_ref(ionstemp,chgtemp,opts)
    ELSE
      chgtemp = chgval
    END IF

    ALLOCATE(dipol(3,ions%nions))
    ALLOCATE(charge_check(ions%nions))
    ALLOCATE(rsquare(ions%nions))
    ALLOCATE(quad(3,3,ions%nions))

    dipol=0.0_q2
    charge_check=0.0_q2
    rsquare=0.0_q2
    quad=0.0_q2

    chgtot = 0.0_q2
    diptot(:) = 0.0_q2
    quadtot(:,:) = 0.0_q2
    rsquaretot = 0.0_q2

    DO n1 = 1,chgval%npts(1)
      DO n2 = 1,chgval%npts(2)
        DO n3 = 1,chgval%npts(3)
         IF (bdr%volnum(n1,n2,n3) == bdr%nvols+1) CYCLE
         tab(1)=REAL(n1,q2)
         tab(2)=REAL(n2,q2)
         tab(3)=REAL(n3,q2)
         ! transform the point to cartesian
         tab_car(:) = lat2car(chgtemp, tab(:))
         ! get the corresponding ion number
         ion_of_point=bdr%nnion(bdr%volnum(n1,n2,n3))
         ! recalculate charge
         charge_check(ion_of_point)=charge_check(ion_of_point)+chgval%rho(n1,n2,n3)
         ! point position
         rr0(:)=tab_car(:)-ions%r_car(ion_of_point,:)
         CALL dpbc_car(ions,rr0)
         distance=SQRT(DOT_PRODUCT(rr0,rr0))
         ! dipole stuff
         dipol(:,ion_of_point)=dipol(:,ion_of_point)+rr0(:)*chgval%rho(n1,n2,n3)
         ! quadrupole stuff
         rsquare(ion_of_point)=rsquare(ion_of_point)+distance*distance*chgval%rho(n1,n2,n3)
         quad(1,1,ion_of_point)=quad(1,1,ion_of_point)+(3.0d0*rr0(1)*rr0(1)-distance*distance)*chgval%rho(n1,n2,n3)
         quad(2,2,ion_of_point)=quad(2,2,ion_of_point)+(3.0d0*rr0(2)*rr0(2)-distance*distance)*chgval%rho(n1,n2,n3)
         quad(3,3,ion_of_point)=quad(3,3,ion_of_point)+(3.0d0*rr0(3)*rr0(3)-distance*distance)*chgval%rho(n1,n2,n3)
         quad(1,2,ion_of_point)=quad(1,2,ion_of_point)+ 3.0d0*rr0(1)*rr0(2)                   *chgval%rho(n1,n2,n3)
         quad(1,3,ion_of_point)=quad(1,3,ion_of_point)+ 3.0d0*rr0(1)*rr0(3)                   *chgval%rho(n1,n2,n3)
         quad(2,3,ion_of_point)=quad(2,3,ion_of_point)+ 3.0d0*rr0(2)*rr0(3)                   *chgval%rho(n1,n2,n3)
        END DO
      END DO
    END DO
  
    dipol = dipol/REAL(chgval%nrho,q2)
    charge_check = charge_check/REAL(chgval%nrho,q2)
    dipol = dipol*charge_sign
    charge_check = charge_check*charge_sign

    quad(2,1,:)=quad(1,2,:)
    quad(3,2,:)=quad(2,3,:)
    quad(3,1,:)=quad(1,3,:)
    quad=quad/REAL(chgval%nrho,q2)
    quad=quad*charge_sign
    rsquare=rsquare/REAL(chgval%nrho,q2)
    rsquare=rsquare*charge_sign

    OPEN(345,FILE='DIP.dat',STATUS='replace',ACTION='write')
    OPEN(346,FILE='QPL.dat',STATUS='replace',ACTION='write')

    IF (opts%in_opt == opts%in_chgcar4 .OR. opts%in_opt == opts%in_chgcar5) THEN
     dipol = dipol/bohr2angs
     quad = quad/(bohr2angs*bohr2angs)
     rsquare = rsquare/(bohr2angs*bohr2angs)
     ions%r_car = ions%r_car/bohr2angs
     WRITE (UNIT=345, FMT=*) 'VASP CALCULATION: UNITS ARE CONVERTED TO atomic units '
     WRITE (UNIT=346, FMT=*) 'VASP CALCULATION: UNITS ARE CONVERTED TO atomic units '
    ELSE
     WRITE (UNIT=345, FMT=*) 'DATA FROM CUBE FILE: all quantities are in atomic units '
     WRITE (UNIT=346, FMT=*) 'DATA FROM CUBE FILE: all quantities are in atomic units '
    ENDIF

    WRITE (UNIT=345, FMT='(A,I4)') ' NIONS =', ions%nions 
    WRITE(345,'(2A)') '----------------------------------------------------------------', &
  &                '-------------------------------------------------------'
    WRITE(345,'(1X,1A,9X,2(1A1,11X),3X,1A1,8X,1A6,4(10X,A4))')  &
  & '#','X','Y','Z','CHARGE','DIPX','DIPY','DIPZ','DTOT'
    WRITE(345,'(2A)') '----------------------------------------------------------------', &
  &                '------------------------------------------------'

    WRITE (UNIT=346, FMT='(A,I4)') ' NIONS =', ions%nions 
    WRITE(346,'(2A)') '----------------------------------------------------------------', &
  &                '----------------------------------------------'
    WRITE(346,'(1X,1A,1X,7(10X,A4))')  &
  & '#','QXX ','QYY ','QZZ ','QXY ','QXZ ','QYZ ','R2  '
    WRITE(346,'(2A)') '----------------------------------------------------------------', &
  &                '---------------------------------------'

    DO n1 = 1,ions%nions
      chgtot = chgtot + charge_check(n1)
      diptot(:) = diptot(:) + dipol(:,n1)     
      quadtot(:,:) = quadtot(:,:) + quad(:,:,n1)     
      rsquaretot = rsquaretot + rsquare(n1)
      WRITE (UNIT=345, FMT='(I3,3F12.6,F14.5,X,4F14.5)') n1,ions%r_car(n1,:),charge_check(n1),dipol(:,n1), &
  &                           SQRT(DOT_PRODUCT(dipol(:,n1),dipol(:,n1)))                                   
      WRITE (UNIT=346, FMT='(I3,X,7F14.5)') n1, quad(1,1,n1),quad(2,2,n1),quad(3,3,n1), &
  &                         quad(1,2,n1),quad(1,3,n1),quad(2,3,n1),rsquare(n1)
    END DO
    WRITE(345,'(2A)') '----------------------------------------------------------------', &
  &                '----------------------------------------------'
    WRITE(346,'(2A)') '----------------------------------------------------------------', &
  &                '---------------------------------------'
    WRITE (UNIT=345, FMT='(A,34X,F14.5,X,4F14.5)') 'TOTAL',chgtot,diptot(:), &
  &                           SQRT(DOT_PRODUCT(diptot,diptot))                                   
    WRITE (UNIT=346, FMT='(A,F13.5,6F14.5)') 'TOTAL', quadtot(1,1),quadtot(2,2),quadtot(3,3), &
  &                         quadtot(1,2),quadtot(1,3),quadtot(2,3),rsquaretot
    WRITE(345,'(2A)') '----------------------------------------------------------------', &
  &                '----------------------------------------------'
    WRITE(346,'(2A)') '----------------------------------------------------------------', &
  &                '---------------------------------------'
    CLOSE(345)
    CLOSE(346)

    DEALLOCATE(dipol)
    DEALLOCATE(charge_check)
    DEALLOCATE(quad)
    DEALLOCATE(rsquare)

  RETURN
  END SUBROUTINE multipole_calc

  END MODULE multipole_mod
