!! Copyright (C) 2018 X. Andrade
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the Lesser GNU General Public License as published by
!! the Free Software Foundation; either version 3, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

module pseudo_oct_m

  implicit none

  private

  public ::                              &
    pseudo_t,                            &
    pseudo_detect_format,                &
    pseudo_init,                         &
    pseudo_end,                          &
    pseudo_type,                         &
    pseudo_format,                       &
    pseudo_valence_charge,               &
    pseudo_mesh_spacing,                 &
    pseudo_mesh_size,                    &
    pseudo_mass,                         &
    pseudo_lmax,                         &
    pseudo_llocal,                       &
    pseudo_nchannels,                    &
    pseudo_nprojectors,                  &
    pseudo_grid,                         &
    pseudo_grid_weights,                 &
    pseudo_local_potential,              &
    pseudo_projector,                    &
    pseudo_has_radial_function,          &
    pseudo_radial_function,              &
    pseudo_radial_potential,             &
    pseudo_has_nlcc,                     &
    pseudo_nlcc_density,                 &
    pseudo_dij,                          &
    pseudo_has_density,                  &
    pseudo_density,                      &
    pseudo_nwavefunctions,               &
    pseudo_wavefunction,                 &
    pseudo_exchange,                     &
    pseudo_correlation,                  &
    pseudo_has_total_angular_momentum,   &
    pseudo_projector_2j,                 &
    pseudo_wavefunction_2j
  
  !the following sets of values have to match with those on base.hpp
  integer, parameter, public ::               &
    PSEUDO_TYPE_ULTRASOFT         = 30,       &
    PSEUDO_TYPE_SEMILOCAL         = 31,       &
    PSEUDO_TYPE_KLEINMAN_BYLANDER = 32,       &
    PSEUDO_TYPE_PAW               = 33

  integer, parameter, public ::                       &
    PSEUDO_STATUS_SUCCESS                      = 0,   &
    PSEUDO_STATUS_FILE_NOT_FOUND               = 455, &
    PSEUDO_STATUS_FORMAT_NOT_SUPPORTED         = 456, &
    PSEUDO_STATUS_UNKNOWN_FORMAT               = 457, &
    PSEUDO_STATUS_UNSUPPORTED_TYPE_ULTRASOFT   = 458, &
    PSEUDO_STATUS_UNSUPPORTED_TYPE_PAW         = 459, &
    PSEUDO_STATUS_UNSUPPORTED_TYPE             = 460

  integer, parameter, public ::                       &
    PSEUDO_FORMAT_FILE_NOT_FOUND             = 773,   &
    PSEUDO_FORMAT_UNKNOWN                    = 774,   &
    PSEUDO_FORMAT_UPF1                       = 775,   &
    PSEUDO_FORMAT_UPF2                       = 776,   &
    PSEUDO_FORMAT_QSO                        = 777,   &
    PSEUDO_FORMAT_PSML                       = 778,   &
    PSEUDO_FORMAT_PSF                        = 779,   &
    PSEUDO_FORMAT_CPI                        = 780,   &
    PSEUDO_FORMAT_FHI                        = 781,   &
    PSEUDO_FORMAT_HGH                        = 782,   &
    PSEUDO_FORMAT_PSP8                       = 783

  ! we only define these values here, the specific functionals are
  ! obtained from libxc
  integer, parameter, public ::                       &
    PSEUDO_EXCHANGE_UNKNOWN                  = -2,    &
    PSEUDO_EXCHANGE_ANY                      = -1

  integer, parameter, public ::                       &
    PSEUDO_CORRELATION_UNKNOWN               = -2,    &
    PSEUDO_CORRELATION_ANY                   = -1
  
  type pseudo_t
    private
    integer(8) :: dummy
  end type pseudo_t

  interface
    
    ! -------------------------------------------------
    
    integer function pseudo_detect_format(filename)
      import :: pseudo_t
      implicit none
      
      character(len=*), intent(in)    :: filename
    end function pseudo_detect_format

    ! -------------------------------------------------
    
    subroutine pseudo_init(pseudo, filename, fmt, ierr)
      import :: pseudo_t
      implicit none
      
      type(pseudo_t),   intent(out)   :: pseudo
      character(len=*), intent(in)    :: filename
      integer,          intent(in)    :: fmt
      integer,          intent(out)   :: ierr
    end subroutine pseudo_init

    ! -------------------------------------------------
    
    subroutine pseudo_end(pseudo)
      import :: pseudo_t
      implicit none

      type(pseudo_t),   intent(inout) :: pseudo
    end subroutine pseudo_end

    ! -------------------------------------------------
    
    integer function pseudo_type(pseudo)
      import :: pseudo_t
      implicit none
      
      type(pseudo_t),   intent(in)    :: pseudo
    end function pseudo_type

    ! -------------------------------------------------
    
    integer function pseudo_format(pseudo)
      import :: pseudo_t
      implicit none
      
      type(pseudo_t),   intent(in)    :: pseudo
    end function pseudo_format

    ! -------------------------------------------------
    
    FLOAT function pseudo_valence_charge(pseudo)
      import :: pseudo_t
      implicit none

      type(pseudo_t),   intent(in)    :: pseudo
    end function pseudo_valence_charge

    ! -------------------------------------------------

    FLOAT function pseudo_mesh_spacing(pseudo)
      import :: pseudo_t
      implicit none
      
      type(pseudo_t),   intent(in)    :: pseudo
    end function pseudo_mesh_spacing

    ! -------------------------------------------------
    
    integer function pseudo_mesh_size(pseudo)
      import :: pseudo_t
      implicit none
      
      type(pseudo_t),   intent(in)    :: pseudo
    end function pseudo_mesh_size

    ! -------------------------------------------------
    
    FLOAT function pseudo_mass(pseudo)
      import :: pseudo_t
      implicit none
      
      type(pseudo_t),   intent(in)    :: pseudo
    end function pseudo_mass

    ! -------------------------------------------------
    
    integer function pseudo_lmax(pseudo)
      import :: pseudo_t
      implicit none
      
      type(pseudo_t),   intent(in)    :: pseudo
    end function pseudo_lmax

    ! -------------------------------------------------
    
    integer function pseudo_llocal(pseudo)
      import :: pseudo_t
      implicit none
      
      type(pseudo_t),   intent(in)    :: pseudo
    end function pseudo_llocal

    ! -------------------------------------------------

    integer function pseudo_nchannels(pseudo)
      import :: pseudo_t
      implicit none
      
      type(pseudo_t),   intent(in)    :: pseudo
    end function pseudo_nchannels

    ! -------------------------------------------------

    integer function pseudo_nprojectors(pseudo)
      import :: pseudo_t
      implicit none
      
      type(pseudo_t),   intent(in)    :: pseudo
    end function pseudo_nprojectors

    ! -------------------------------------------------

    subroutine pseudo_grid(pseudo, grid)
      import :: pseudo_t
      implicit none
      
      type(pseudo_t),   intent(in)    :: pseudo
      FLOAT,            intent(out)   :: grid
    end subroutine pseudo_grid
    
    ! -------------------------------------------------

    subroutine pseudo_grid_weights(pseudo, weight)
      import :: pseudo_t
      implicit none
      
      type(pseudo_t),   intent(in)    :: pseudo
      FLOAT,            intent(out)   :: weight
    end subroutine pseudo_grid_weights

    ! -------------------------------------------------

    subroutine pseudo_local_potential(pseudo, local_potential)
      import :: pseudo_t
      implicit none
      
      type(pseudo_t),   intent(in)    :: pseudo
      FLOAT,            intent(out)   :: local_potential
    end subroutine pseudo_local_potential
    
    ! -------------------------------------------------
    
    subroutine pseudo_projector(pseudo, l, ic, projector)
      import :: pseudo_t
      implicit none
      
      type(pseudo_t),   intent(in)    :: pseudo
      integer,          intent(in)    :: l
      integer,          intent(in)    :: ic
      FLOAT,            intent(out)   :: projector
    end subroutine pseudo_projector

    ! -------------------------------------------------
    
    FLOAT function pseudo_dij(pseudo, l, ic, jc)
      import :: pseudo_t
      implicit none
      
      type(pseudo_t),   intent(in)    :: pseudo
      integer,          intent(in)    :: l
      integer,          intent(in)    :: ic
      integer,          intent(in)    :: jc
    end function pseudo_dij

    ! -------------------------------------------------
    
    subroutine pseudo_radial_potential(pseudo, l, radial_potential)
      import :: pseudo_t
      implicit none
      
      type(pseudo_t),   intent(in)    :: pseudo
      integer,          intent(in)    :: l
      FLOAT,            intent(inout) :: radial_potential
    end subroutine pseudo_radial_potential

    ! -------------------------------------------------
    
    subroutine pseudo_radial_function(pseudo, l, radial_function)
      import :: pseudo_t
      implicit none
      
      type(pseudo_t),   intent(in)    :: pseudo
      integer,          intent(in)    :: l
      FLOAT,            intent(out)   :: radial_function
    end subroutine pseudo_radial_function

    ! -------------------------------------------------
    
    subroutine pseudo_nlcc_density(pseudo, nlcc_density)
      import :: pseudo_t
      implicit none
      
      type(pseudo_t),   intent(in)    :: pseudo
      FLOAT,            intent(out)   :: nlcc_density
    end subroutine pseudo_nlcc_density

    ! -------------------------------------------------
    
    subroutine pseudo_density(pseudo, density)
      import :: pseudo_t
      implicit none
      
      type(pseudo_t),   intent(in)    :: pseudo
      FLOAT,            intent(out)   :: density
    end subroutine pseudo_density

    ! -------------------------------------------------
    
    integer function pseudo_nwavefunctions(pseudo)
      import :: pseudo_t
      implicit none
      
      type(pseudo_t),   intent(in)    :: pseudo
    end function pseudo_nwavefunctions

    ! -------------------------------------------------
    
    subroutine pseudo_wavefunction(pseudo, index, n, l, occ, wf)
      import :: pseudo_t
      implicit none
      
      type(pseudo_t),   intent(in)    :: pseudo
      integer,          intent(in)    :: index
      integer,          intent(out)   :: n
      integer,          intent(out)   :: l
      FLOAT,            intent(out)   :: occ
      FLOAT,            intent(out)   :: wf
    end subroutine pseudo_wavefunction

    ! -------------------------------------------------

    integer function pseudo_exchange(pseudo)
      import :: pseudo_t
      implicit none
      
      type(pseudo_t),   intent(in)    :: pseudo
    end function pseudo_exchange

    ! -------------------------------------------------

    integer function pseudo_correlation(pseudo)
      import :: pseudo_t
      implicit none
      
      type(pseudo_t),   intent(in)    :: pseudo
    end function pseudo_correlation

    ! -------------------------------------------------
    
    integer function pseudo_projector_2j(pseudo, l, ic)
      import :: pseudo_t
      implicit none
      
      type(pseudo_t),   intent(in)    :: pseudo
      integer,          intent(in)    :: l
      integer,          intent(in)    :: ic
    end function pseudo_projector_2j

    ! -------------------------------------------------
    
    integer function pseudo_wavefunction_2j(pseudo, ii)
      import :: pseudo_t
      implicit none
      
      type(pseudo_t),   intent(in)    :: pseudo
      integer,          intent(in)    :: ii
    end function pseudo_wavefunction_2j
    
  end interface
  
contains

  ! -------------------------------------------------
  
  logical function pseudo_has_nlcc(pseudo)
    type(pseudo_t),   intent(in)      :: pseudo
    
    interface
      integer function pseudo_has_nlcc_low(pseudo)
        import :: pseudo_t
        implicit none
        
        type(pseudo_t),   intent(in)      :: pseudo
      end function pseudo_has_nlcc_low
    end interface

    pseudo_has_nlcc = (pseudo_has_nlcc_low(pseudo) /= 0)
    
  end function pseudo_has_nlcc

  ! -------------------------------------------------
  
  logical function pseudo_has_density(pseudo)
    type(pseudo_t),   intent(in)      :: pseudo
    
    interface
      integer function pseudo_has_density_low(pseudo)
        import :: pseudo_t
        implicit none
        
        type(pseudo_t),   intent(in)      :: pseudo
      end function pseudo_has_density_low
    end interface

    pseudo_has_density = (pseudo_has_density_low(pseudo) /= 0)
    
  end function pseudo_has_density
    ! -------------------------------------------------
  
  logical function pseudo_has_total_angular_momentum(pseudo)
    type(pseudo_t),   intent(in)      :: pseudo
    
    interface
      integer function pseudo_has_total_angular_momentum_low(pseudo)
        import :: pseudo_t
        implicit none
        
        type(pseudo_t),   intent(in)      :: pseudo
      end function pseudo_has_total_angular_momentum_low
    end interface

    pseudo_has_total_angular_momentum = (pseudo_has_total_angular_momentum_low(pseudo) /= 0)
    
  end function pseudo_has_total_angular_momentum

  ! -------------------------------------------------
  
  logical function pseudo_has_radial_function(pseudo, l)
    type(pseudo_t),   intent(in)      :: pseudo
    integer,          intent(in)      :: l
    
    interface
      integer function pseudo_has_radial_function_low(pseudo, l)
        import :: pseudo_t
        implicit none
        
        type(pseudo_t),   intent(in)      :: pseudo
        integer,          intent(in)      :: l
      end function pseudo_has_radial_function_low
    end interface
    
    pseudo_has_radial_function = (pseudo_has_radial_function_low(pseudo, l) /= 0)
    
  end function pseudo_has_radial_function

end module pseudo_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
