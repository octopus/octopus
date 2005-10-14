!! Copyright (C) 2002-2003 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
!! $Id$

#include "global.h"

module mesh_function
  use global
  use messages
  use lib_basic_alg
  use math
  use mesh
  use cube_function
  use par_vec
  use profiling_mod
#ifdef HAVE_MPI
  use mpi_mod
#endif
#if defined(HAVE_MPI) && !defined(MPI_H)
  use mpi
#endif


  implicit none


#if defined(HAVE_MPI) && defined(MPI_H)
# include "mpif.h"
#endif

  private
  public :: dmf_integrate, zmf_integrate, &
            dmf_dotp, zmf_dotp, &
            dmf_nrm2, zmf_nrm2, &
            dmf_moment, zmf_moment, &
            dmf_random, zmf_random


contains

#include "undef.F90"
#include "real.F90"
#include "mf_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "mf_inc.F90"

end module mesh_function
