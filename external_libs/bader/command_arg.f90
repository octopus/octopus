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

! conversion between
!    COMMAND_ARGUMENT_COUNT and IARGC
!    GET_COMMAND_ARGUMENT and GETARG

  FUNCTION COMMAND_ARGUMENT_COUNT()
    INTEGER :: COMMAND_ARGUMENT_COUNT
    COMMAND_ARGUMENT_COUNT=IARGC()
    RETURN
  END FUNCTION

  SUBROUTINE GET_COMMAND_ARGUMENT(m,p)
    CHARACTER(LEN=128) :: p
    INTEGER :: m
    CALL GETARG(m,p)
  END SUBROUTINE
