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
!  Module with ion data structure
!-----------------------------------------------------------------------------------!

MODULE ions_mod
  USE kind_mod
  IMPLICIT NONE

  TYPE :: ions_obj
    REAL(q2),ALLOCATABLE,DIMENSION(:,:) :: r_car,r_dir,r_lat
    REAL(q2),ALLOCATABLE,DIMENSION(:) :: ion_chg
    REAL(q2),DIMENSION(3,3) :: lattice,dir2car,car2dir
    INTEGER,ALLOCATABLE,DIMENSION(:) :: num_ion
    INTEGER,ALLOCATABLE,DIMENSION(:) :: atomic_num
    CHARACTER*330:: name_ion
    INTEGER :: niontypes,nions
  END TYPE

  PRIVATE
  PUBLIC :: ions_obj

!-----------------------------------------------------------------------------------!

END MODULE ions_mod
