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

!-----------------------------------------------------------------------------------!
! Bader charge density analysis program
!  Module for analyzing the charge
!-----------------------------------------------------------------------------------!

MODULE charge_mod

  USE kind_mod
  USE matrix_mod
  USE ions_mod
  IMPLICIT NONE

  !Q define weight type.
  TYPE, PUBLIC :: weight
    real(q2),ALLOCATABLE, DIMENSION(:) :: w
  END TYPE
  !Q
  TYPE :: charge_obj
    TYPE(weight), ALLOCATABLE, DIMENSION(:,:,:) :: weight !Q
    REAL(q2),ALLOCATABLE,DIMENSION(:,:,:) :: rho
    REAL(q2),DIMENSION(3,3) :: lat2car,car2lat
    REAL(q2),DIMENSION(-1:1,-1:1,-1:1) :: lat_dist, lat_i_dist
    REAL(q2),DIMENSION(3) :: org_lat, org_car !,org_dir
    REAL(q2),DIMENSION(3) :: i_npts
    INTEGER,DIMENSION(3) :: npts
    INTEGER :: nrho
  END TYPE

  PRIVATE
  PUBLIC :: charge_obj
  PUBLIC :: rho_val, rho_grad, rho_grad_dir
  PUBLIC :: pbc, dpbc_dir, dpbc, dpbc_dir_org, dpbc_car, pbc_r_lat
  PUBLIC :: to_lat, is_max, is_max_ongrid
  PUBLIC :: lat2car, car2lat, lat2dir, dir2lat

  INTERFACE ASSIGNMENT(=)
    MODULE PROCEDURE assign_charge
  END INTERFACE

  CONTAINS

!-----------------------------------------------------------------------------------!
!  assign: copy one charge object to another
!-----------------------------------------------------------------------------------!

  SUBROUTINE assign_charge(chg1,chg2)

    TYPE(charge_obj), INTENT(INOUT) :: chg1
    TYPE(charge_obj), INTENT(IN) :: chg2
    
    chg1%lat2car = chg2%lat2car
    chg1%car2lat = chg2%car2lat
    chg1%org_lat = chg2%org_lat
!    chg1%org_dir = chg2%org_dir
    chg1%org_car = chg2%org_car
    chg1%lat_dist = chg2%lat_dist
    chg1%lat_i_dist = chg2%lat_i_dist
    chg1%i_npts = chg2%i_npts
    chg1%npts = chg2%npts
    chg1%nrho = chg2%nrho

    ALLOCATE(chg1%rho(chg1%npts(1),chg1%npts(2),chg1%npts(3)))
    chg1%rho=chg2%rho

    END SUBROUTINE

!-----------------------------------------------------------------------------------!
!  rho_val:  Return the density at the point (p1,p2,p3) taking into account the
!    boundary conditions.  This function is used to address points outside the
!    charge density array without a bunch of if statements at the place the value
!    is needed.
!-----------------------------------------------------------------------------------!

  FUNCTION rho_val(chg,p1,p2,p3)

    TYPE(charge_obj) :: chg
    INTEGER,INTENT(IN) :: p1, p2, p3
    REAL(q2) :: rho_val

    INTEGER,DIMENSION(3) :: p
    INTEGER :: i

    p=(/p1,p2,p3/)
    DO i=1,3
      DO
        IF(p(i) >= 1) EXIT
        p(i) = p(i) + chg%npts(i)
      END DO
      DO
        IF(p(i) <= chg%npts(i)) EXIT
        p(i) = p(i) - chg%npts(i)
      END DO
    END DO

    rho_val=chg%rho(p(1),p(2),p(3))

  RETURN
  END FUNCTION rho_val

!-----------------------------------------------------------------------------------!
!  rho_grad:  Return the density and gradient at the point r
!-----------------------------------------------------------------------------------!

  FUNCTION rho_grad(chg,r,rho)

    TYPE(charge_obj) :: chg
    REAL(q2),DIMENSION(3),INTENT(IN) :: r
    REAL(q2),INTENT(OUT) :: rho
    REAL(q2),DIMENSION(3) :: rho_grad

    INTEGER :: p1, p2, p3
    REAL(q2),DIMENSION(3) :: rho_grad_lat
    REAL(q2) :: f1, f2, f3, g1, g2, g3
    REAL(q2) :: rho000, rho001, rho010, rho100, rho011, rho101, rho110, rho111
    REAL(q2) :: rho00_, rho01_, rho10_, rho11_
    REAL(q2) :: rho0__, rho1__, rho_0_, rho_1_, rho__0, rho__1
    REAL(q2) :: rho_00, rho_01, rho_10, rho_11

    p1 = FLOOR(r(1))
    p2 = FLOOR(r(2))
    p3 = FLOOR(r(3))

    f1 = r(1) - REAL(p1,q2)
    f2 = r(2) - REAL(p2,q2)
    f3 = r(3) - REAL(p3,q2)

    g1 = 1._q2-f1
    g2 = 1._q2-f2
    g3 = 1._q2-f3

    rho000 = rho_val(chg,p1,p2,p3)
    rho001 = rho_val(chg,p1,p2,p3+1)
    rho010 = rho_val(chg,p1,p2+1,p3)
    rho100 = rho_val(chg,p1+1,p2,p3)
    rho011 = rho_val(chg,p1,p2+1,p3+1)
    rho101 = rho_val(chg,p1+1,p2,p3+1)
    rho110 = rho_val(chg,p1+1,p2+1,p3)
    rho111 = rho_val(chg,p1+1,p2+1,p3+1)

!    write(*,'(A,2F12.4)')'   011, 111: ',rho011,rho111
!    write(*,'(A,2F12.4)')'   001, 101: ',rho001,rho101
!    write(*,'(A,2F12.4)')'   010, 110: ',rho010,rho110
!    write(*,'(A,2F12.4)')'   000, 100: ',rho000,rho100

    rho00_ = rho000*g3 + rho001*f3
    rho01_ = rho010*g3 + rho011*f3
    rho10_ = rho100*g3 + rho101*f3
    rho11_ = rho110*g3 + rho111*f3

    rho0__ = rho00_*g2 + rho01_*f2
    rho1__ = rho10_*g2 + rho11_*f2

    rho = rho0__*g1 + rho1__*f1

! More work for gradients

    rho_0_ = rho00_*g1 + rho10_*f1
    rho_1_ = rho01_*g1 + rho11_*f1

    rho_00 = rho000*g1 + rho100*f1
    rho_01 = rho001*g1 + rho101*f1
    rho_10 = rho010*g1 + rho110*f1
    rho_11 = rho011*g1 + rho111*f1

    rho__0 = rho_00*g2 + rho_10*f2
    rho__1 = rho_01*g2 + rho_11*f2

    rho_grad_lat(1) = rho1__ - rho0__
    rho_grad_lat(2) = rho_1_ - rho_0_
    rho_grad_lat(3) = rho__1 - rho__0

    CALL vector_matrix(rho_grad_lat, chg%car2lat, rho_grad)

  RETURN
  END FUNCTION rho_grad

!-----------------------------------------------------------------------------------!
!  rho_grad_dir:  Return the direction of the gradient in lattice vectors
!                 at the grid position p
!-----------------------------------------------------------------------------------!

  FUNCTION rho_grad_dir(chg,p)

    TYPE(charge_obj) :: chg
    INTEGER,DIMENSION(3),INTENT(IN) :: p
    REAL(q2),DIMENSION(3) :: rho_grad_dir

    INTEGER :: p1, p2, p3
    REAL(q2),DIMENSION(3) :: rho_grad_lat, rho_grad_car
    REAL(q2) :: rho000, rho001, rho010, rho100, rho00_1, rho_100, rho0_10

    p1 = p(1)
    p2 = p(2)
    p3 = p(3)
    
    rho000 = rho_val(chg,p1,p2,p3)
    rho001 = rho_val(chg,p1,p2,p3+1) 
    rho010 = rho_val(chg,p1,p2+1,p3)
    rho100 = rho_val(chg,p1+1,p2,p3)
    rho00_1 = rho_val(chg,p1,p2,p3-1)
    rho_100 = rho_val(chg,p1-1,p2,p3)
    rho0_10 = rho_val(chg,p1,p2-1,p3)

    rho_grad_lat(1) = (rho100-rho_100)/2._q2
    rho_grad_lat(2) = (rho010-rho0_10)/2._q2
    rho_grad_lat(3) = (rho001-rho00_1)/2._q2

    IF(rho100 < rho000.AND.rho_100 < rho000) rho_grad_lat(1) = 0._q2
    IF(rho010 < rho000.AND.rho0_10 < rho000) rho_grad_lat(2) = 0._q2
    IF(rho001 < rho000.AND.rho00_1 < rho000) rho_grad_lat(3) = 0._q2

    ! convert to cartesian coordinates
    CALL vector_matrix(rho_grad_lat,chg%car2lat,rho_grad_car)
!    CALL matrix_vector(chg%lat2car,rho_grad_lat,rho_grad_car)

    ! express this vector in direct coordinates
    CALL matrix_vector(chg%car2lat,rho_grad_car,rho_grad_dir)

    ! return a unit vector
!    rho_grad_dir=rho_grad_dir/SQRT(SUM(rho_grad_dir*rho_grad_dir))

  RETURN
  END FUNCTION rho_grad_dir
!-----------------------------------------------------------------------------------!
! pbc: Wrap the point (p(1),p(2),p(3)) to the boundary conditions [0,pmax].
!-----------------------------------------------------------------------------------!

  SUBROUTINE pbc_r_lat(r_lat,pmax)

    REAL(q2),DIMENSION(3),INTENT(INOUT) :: r_lat
    INTEGER,DIMENSION(3),INTENT(IN) :: pmax

    INTEGER :: i

    DO i=1,3
      DO
        IF(r_lat(i) > 0) EXIT
        r_lat(i)=r_lat(i)+pmax(i)
      END DO
      DO
        IF(r_lat(i) <= REAL(pmax(i))) EXIT
        r_lat(i)=r_lat(i)-pmax(i)
      END DO
    END DO

  RETURN
  END SUBROUTINE pbc_r_lat

!-----------------------------------------------------------------------------------!
! pbc: Wrap the point (p(1),p(2),p(3)) to the boundary conditions [0,pmax].
!-----------------------------------------------------------------------------------!

  SUBROUTINE pbc(p,pmax)

    INTEGER,DIMENSION(3),INTENT(INOUT) :: p
    INTEGER,DIMENSION(3),INTENT(IN) :: pmax

    INTEGER :: i

    DO i=1,3
      DO
        IF(p(i) > 0) EXIT
        p(i) = p(i) + pmax(i)
      END DO
      DO
        IF(p(i) <= pmax(i)) EXIT
        p(i) = p(i) - pmax(i)
      END DO
    END DO

  RETURN
  END SUBROUTINE pbc

!-----------------------------------------------------------------------------------!
! dpbc:  Wrap the vector dr to the boundary conditions [-ngf/2,ngf/2].
!-----------------------------------------------------------------------------------!

  SUBROUTINE dpbc(dr,nf,nf_2)

    REAL(q2),INTENT(IN),DIMENSION(3) :: nf, nf_2
    REAL(q2),INTENT(INOUT),DIMENSION(3) :: dr

    INTEGER :: i

    DO i=1,3
      DO
        IF(dr(i) > -nf_2(i)) EXIT
        dr(i) = dr(i) + nf(i)
      END DO
      DO
        IF(dr(i) < nf_2(i)) EXIT
        dr(i) = dr(i) - nf(i)
      END DO
    END DO

  RETURN
  END SUBROUTINE dpbc
!-----------------------------------------------------------------------------------!
! dpbc_dir_org:  Wrap the vector dr to the boundary conditions [-1/2,1/2].
!-----------------------------------------------------------------------------------!

  SUBROUTINE dpbc_dir_org(dr)

    REAL(q2),INTENT(INOUT),DIMENSION(3) :: dr

    INTEGER :: i

    DO i=1,3
      DO
        IF(dr(i) > -0.5_q2) EXIT
        dr(i) = dr(i) + 1.0_q2
      END DO
      DO
        IF(dr(i) < 0.5_q2) EXIT
        dr(i) = dr(i) - 1.0_q2
      END DO
    END DO
  RETURN
  END SUBROUTINE dpbc_dir_org

!-----------------------------------------------------------------------------------!
! dpbc_dir:  Find the smallest distance vector in direct coordinates
!-----------------------------------------------------------------------------------!

  SUBROUTINE dpbc_dir(ions,dr_dir)

    TYPE(ions_obj) :: ions
    REAL(q2),INTENT(INOUT),DIMENSION(3) :: dr_dir
    REAL(q2),DIMENSION(3) :: dr_car

    ! convert to cartesian coordinates
    CALL matrix_vector(ions%dir2car, dr_dir, dr_car)

    CALL dpbc_car(ions, dr_car)

    ! express this vector in direct coordinates
    CALL matrix_vector(ions%car2dir, dr_car, dr_dir)

  RETURN
  END SUBROUTINE dpbc_dir

!-----------------------------------------------------------------------------------!
! dpbc_car:  Find the smallest distance vector in Cartesian coordinates
!-----------------------------------------------------------------------------------!

  SUBROUTINE dpbc_car(ions,dr_car)

    TYPE(ions_obj) :: ions
    REAL(q2),INTENT(INOUT),DIMENSION(3) :: dr_car
    REAL(q2),DIMENSION(3) :: drt_car,v1,v2,v3
    REAL(q2) :: dsq,dsqmin
    INTEGER :: d1,d2,d3
    LOGICAL :: done

    dsqmin=DOT_PRODUCT(dr_car,dr_car)
    DO
      done = .TRUE.
      DO d1=-1,1
        v1 = ions%lattice(1,:)*REAL(d1,q2)
        DO d2=-1,1
          v2 = ions%lattice(2,:)*REAL(d2,q2)
          DO d3=-1,1
            v3 = ions%lattice(3,:)*REAL(d3,q2)

            drt_car = dr_car+v1+v2+v3
            dsq = DOT_PRODUCT(drt_car,drt_car)
            IF(dsq<dsqmin) THEN
              dr_car = drt_car
              dsqmin = dsq
              done = .FALSE.
            END IF

          END DO
        END DO
      END DO
      IF(done) EXIT
    END DO

  RETURN
  END SUBROUTINE dpbc_car

!-----------------------------------------------------------------------------------!
! to_lat: return the nearest (integer) lattice  point p to the (read) point r
!-----------------------------------------------------------------------------------!

  FUNCTION to_lat(chg,r)

    TYPE(charge_obj) :: chg
    REAL(q2),DIMENSION(3),INTENT(IN) :: r
    INTEGER,DIMENSION(3) :: to_lat

    INTEGER,DIMENSION(3) :: p0, p, pmin
    REAL(q2),DIMENSION(3) :: f, d_lat, d_car
    REAL(q2) :: dsq, dsq_min
    INTEGER p1, p2, p3
    LOGICAL :: init_flag = .TRUE.

    p0 = FLOOR(r)
    pmin = (/0,0,0/)
    f = r - REAL(p0,q2)
    dsq_min = 0._q2

!    write(*,*) 'into to_lat'
    DO p1=0,1
      DO p2=0,1
        DO p3=0,1
          p=(/p1,p2,p3/)
!          write(*,*) 'testing p',p
          d_lat = REAL(p,q2)-f
          CALL matrix_vector(chg%lat2car, d_lat, d_car)
          dsq = SUM(d_car*d_car)
          IF ((dsq<dsq_min) .OR. init_flag) THEN
            init_flag = .FALSE.
            pmin = p
!            write(*,'(A,3I4,F12.4)') '   pmin: ',p1,p2,p3,dsq
            dsq_min = dsq
          END IF
        END DO
      END DO
    END DO
    to_lat = p0 + pmin
    CALL pbc(to_lat, chg%npts)
!    write(*,*) ' pmin: ',pmin

!    write(*,*) 'out to_lat'

  RETURN
  END FUNCTION to_lat

!-----------------------------------------------------------------------------------!
! is_max: return .true. if the grid point is a maximum of charge density
!-----------------------------------------------------------------------------------!

  FUNCTION is_max(chg,p)

    TYPE(charge_obj) :: chg
    INTEGER,DIMENSION(3),INTENT(IN) :: p
    LOGICAL :: is_max

    REAL(q2) :: rho
    INTEGER :: d1, d2, d3, p1, p2, p3
  
    is_max=.TRUE. 
    p1 = p(1)
    p2 = p(2)
    p3 = p(3)
    rho=rho_val(chg,p1,p2,p3)
!    write(*,*) 'rho: ',rho
    DO d1=-1,1
      p1=p(1)+d1
      DO d2=-1,1
        p2=p(2)+d2
        DO d3=-1,1
          p3=p(3)+d3
!          write(*,*) ' ',d1,d2,d3,rho_val(chg,p1,p2,p3)
          IF(rho_val(chg,p1,p2,p3) > rho) THEN
!            write(*,*) ' false'
            is_max = .FALSE.
          END IF
        END DO
      END DO
    END DO

  RETURN 
  END FUNCTION is_max

!-----------------------------------------------------------------------------------!
! is_max_ongrid: return .true. if the grid point is a maximum of charge density
!-----------------------------------------------------------------------------------!

  FUNCTION is_max_ongrid(chg,p)

    TYPE(charge_obj) :: chg
    INTEGER,DIMENSION(3),INTENT(IN) :: p
    LOGICAL :: is_max_ongrid
    INTEGER :: p1, p2, p3, i
    REAL(q2) :: rho000, rho001, rho010, rho100, rho00_1, rho_100, rho0_10

    p1 = p(1)
    p2 = p(2)
    p3 = p(3)
    is_max_ongrid = .FALSE.

! see which of these is more efficient:
    rho000=rho_val(chg,p1,p2,p3)
    IF(rho_val(chg,p1,p2,p3+1) > rho000) RETURN
    IF(rho_val(chg,p1,p2,p3-1) > rho000) RETURN
    IF(rho_val(chg,p1,p2+1,p3) > rho000) RETURN
    IF(rho_val(chg,p1,p2-1,p3) > rho000) RETURN
    IF(rho_val(chg,p1+1,p2,p3) > rho000) RETURN
    IF(rho_val(chg,p1-1,p2,p3) > rho000) RETURN
    is_max_ongrid = .TRUE.


  RETURN
  END FUNCTION is_max_ongrid

!-----------------------------------------------------------------------------------!
! lat2dir: convert from a lattice grid point to a direct coordinate vector
!-----------------------------------------------------------------------------------!

  FUNCTION lat2dir(chg,p)

    TYPE(charge_obj) :: chg
    REAL(q2),DIMENSION(3),INTENT(IN) :: p
    REAL(q2),DIMENSION(3) :: lat2dir

    lat2dir = p - chg%org_lat
    lat2dir = lat2dir*chg%i_npts
!    lat2dir=lat2dir+chg%org_dir

  RETURN
  END FUNCTION lat2dir

!-----------------------------------------------------------------------------------!
! lat2car: convert from a lattice grid point to a Cartesian coordinate vector
!-----------------------------------------------------------------------------------!

  FUNCTION lat2car(chg,p)

    TYPE(charge_obj) :: chg
    REAL(q2),DIMENSION(3),INTENT(IN) :: p
    REAL(q2),DIMENSION(3) :: lat2car, v

    v = p - chg%org_lat
    CALL matrix_vector(chg%lat2car,v,lat2car)
    lat2car = lat2car + chg%org_car

  RETURN
  END FUNCTION lat2car

!-----------------------------------------------------------------------------------!
! dir2lat: convert from a direct coordinate vector to a lattice grid point
!-----------------------------------------------------------------------------------!

  FUNCTION dir2lat(chg,p)

    TYPE(charge_obj) :: chg
    REAL(q2),DIMENSION(3),INTENT(IN) :: p
    REAL(q2),DIMENSION(3) :: dir2lat

!    dir2lat=p-chg%org_dir
    dir2lat = dir2lat*REAL(chg%npts,q2)
!GH error:    dir2lat=dir2lat+chg%org_dir
    dir2lat = dir2lat + chg%org_lat

  RETURN
  END FUNCTION dir2lat

!-----------------------------------------------------------------------------------!
! car2lat: convert from a Cartesian coordinate vector to a lattice grid point
!-----------------------------------------------------------------------------------!

  FUNCTION car2lat(chg,p)

    TYPE(charge_obj) :: chg
    REAL(q2),DIMENSION(3),INTENT(IN) :: p
    REAL(q2),DIMENSION(3) :: car2lat, v

    v = p - chg%org_car
    CALL matrix_vector(chg%car2lat, v, car2lat)
    car2lat = car2lat + chg%org_lat

  RETURN
  END FUNCTION car2lat

!-----------------------------------------------------------------------------------!

END MODULE charge_mod
