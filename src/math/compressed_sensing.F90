!! Copyright (C) 2011 X. Andrade
!!
!! This program is free software: you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation, either version 3 of the License, or
!! (at your option) any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!!
!! $Id: compressed_sensing.F90 7842 2011-05-05 16:20:39Z joseba $
  
#include <global.h>

module compressed_sensing_m
  use bpdn_m
  use global_m
  use messages_m
  use profiling_m

  implicit none

  private
  
  public ::                                        &
    compressed_sensing_t,                          &
    compressed_sensing_init,                       &
    compressed_sensing_end,                        &
    compressed_sensing_spectral_analysis

  type compressed_sensing_t
    FLOAT   :: sigma
    integer :: ntime
    FLOAT   :: dtime
    FLOAT   :: stime
    integer :: nfreq
    FLOAT   :: dfreq
    FLOAT   :: sfreq
    FLOAT, pointer :: fourier_matrix(:, :)
  end type compressed_sensing_t

contains

  subroutine compressed_sensing_init(this, ntime, dtime, stime, nfreq, dfreq, sfreq)
    type(compressed_sensing_t),  intent(out) :: this
    integer,                     intent(in)  :: ntime
    FLOAT,                       intent(in)  :: dtime
    FLOAT,                       intent(in)  :: stime
    integer,                     intent(in)  :: nfreq
    FLOAT,                       intent(in)  :: dfreq
    FLOAT,                       intent(in)  :: sfreq

    integer :: itime, ifreq
    FLOAT   :: time, freq

    PUSH_SUB(compressed_sensing_init)

    this%sigma = CNST(3.0e-4)

    this%ntime = ntime
    this%dtime = dtime
    this%stime = stime
    this%nfreq = nfreq
    this%dfreq = dfreq
    this%sfreq = sfreq
    
    SAFE_ALLOCATE(this%fourier_matrix(1:this%ntime, 1:this%nfreq))

    do ifreq = 1, this%nfreq
      freq = ifreq*this%dfreq + this%sfreq
      do itime = 1, this%ntime
        time = itime*this%dtime + this%stime

        this%fourier_matrix(itime, ifreq) = sin(freq*time)
      end do
    end do
    
    POP_SUB(compressed_sensing_init)
  end subroutine compressed_sensing_init

  ! -------------------------------------------------------------------

  subroutine compressed_sensing_end(this)
    type(compressed_sensing_t),  intent(inout) :: this

    PUSH_SUB(compressed_sensing_end)

    SAFE_DEALLOCATE_P(this%fourier_matrix)

    POP_SUB(compressed_sensing_end)
  end subroutine compressed_sensing_end
  
  ! -------------------------------------------------------------------

  subroutine compressed_sensing_spectral_analysis(this, time_function, freq_function)
    type(compressed_sensing_t),  intent(out) :: this
    FLOAT,                       intent(in)  :: time_function(:)
    FLOAT,                       intent(out) :: freq_function(:)
    
    PUSH_SUB(compressed_sensing_spectral_analysis)

    call bpdn(this%ntime, this%nfreq, this%fourier_matrix, time_function, this%sigma, freq_function)

    POP_SUB(compressed_sensing_spectral_analysis)    
  end subroutine compressed_sensing_spectral_analysis

end module compressed_sensing_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
