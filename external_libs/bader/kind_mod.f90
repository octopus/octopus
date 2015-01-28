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
!   Module specifying precision and constants
!-----------------------------------------------------------------------------------!
MODULE kind_mod
  IMPLICIT NONE

  PUBLIC

! Public parameters
  INTEGER,PARAMETER :: q1=SELECTED_REAL_KIND(6,30)         ! single precision
  INTEGER,PARAMETER :: q2=SELECTED_REAL_KIND(15,305)       ! double precision
  REAL(q2),PARAMETER :: pi=3.141592653589793238462643_q2    

END MODULE kind_mod

