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
  use blas_m
  use bpdn_m
  use global_m
  use messages_m
  use profiling_m

  implicit none

  private
  
  ! these values are copied from td/spectrum.f90
  integer, parameter ::            &
    SPECTRUM_TRANSFORM_EXP   = 1,  &
    SPECTRUM_TRANSFORM_SIN   = 2,  &
    SPECTRUM_TRANSFORM_COS   = 3
  
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
    type(bpdn_matrix) :: fourier_matrix
  end type compressed_sensing_t

contains

  subroutine compressed_sensing_init(this, transform_type, ntime, dtime, stime, nfreq, dfreq, sfreq, noise)
    type(compressed_sensing_t),  intent(out) :: this
    integer,                     intent(in)  :: transform_type
    integer,                     intent(in)  :: ntime
    FLOAT,                       intent(in)  :: dtime
    FLOAT,                       intent(in)  :: stime
    integer,                     intent(in)  :: nfreq
    FLOAT,                       intent(in)  :: dfreq
    FLOAT,                       intent(in)  :: sfreq
    FLOAT,                       intent(in)  :: noise

    integer :: itime, ifreq, type
    FLOAT   :: time, freq

    PUSH_SUB(compressed_sensing_init)

    this%sigma = noise

    this%ntime = ntime
    this%dtime = dtime
    this%stime = stime
    this%nfreq = nfreq
    this%dfreq = dfreq
    this%sfreq = sfreq
    
    if(transform_type == SPECTRUM_TRANSFORM_EXP) then
      
      call bpdn_matrix_init(this%fourier_matrix, this%ntime, this%nfreq, EXPLICIT_MATRIX)
      
      do ifreq = 1, this%nfreq
        freq = (ifreq - 1)*this%dfreq + this%sfreq
        
        select case(transform_type)
        case(SPECTRUM_TRANSFORM_EXP)
          do itime = 1, this%ntime
            time = (itime - 1)*this%dtime + this%stime
            this%fourier_matrix%matrix(itime, ifreq) = exp(-freq*time)
          end do
          
        case(SPECTRUM_TRANSFORM_SIN)
          do itime = 1, this%ntime
            time = (itime - 1)*this%dtime + this%stime
            this%fourier_matrix%matrix(itime, ifreq) = sin(freq*time)
          end do
          
        case(SPECTRUM_TRANSFORM_COS)
          do itime = 1, this%ntime
            time = (itime - 1)*this%dtime + this%stime
            this%fourier_matrix%matrix(itime, ifreq) = cos(freq*time)
          end do
        end select
        
      end do

    else

      select case(transform_type)
      case(SPECTRUM_TRANSFORM_SIN)
        type = SIN_MATRIX
      case(SPECTRUM_TRANSFORM_COS)
        type = COS_MATRIX
      end select

    call bpdn_matrix_init(this%fourier_matrix, this%ntime, this%nfreq, type)
    call bpdn_matrix_set_delta(this%fourier_matrix, this%dtime, this%dfreq)

  endif
    
    POP_SUB(compressed_sensing_init)
  end subroutine compressed_sensing_init

  ! -------------------------------------------------------------------

  subroutine compressed_sensing_end(this)
    type(compressed_sensing_t),  intent(inout) :: this

    PUSH_SUB(compressed_sensing_end)

    call bpdn_matrix_end(this%fourier_matrix)

    POP_SUB(compressed_sensing_end)
  end subroutine compressed_sensing_end
  
  ! -------------------------------------------------------------------

  subroutine compressed_sensing_spectral_analysis(this, time_function, freq_function)
    type(compressed_sensing_t),  intent(out) :: this
    FLOAT,                       intent(in)  :: time_function(:)
    FLOAT,                       intent(out) :: freq_function(:)

    integer :: ierr
    FLOAT, allocatable :: tf_normalized(:)
    FLOAT :: nrm
    
    PUSH_SUB(compressed_sensing_spectral_analysis)

    SAFE_ALLOCATE(tf_normalized(1:this%ntime))

    ! to avoid numerical problems we work with a normalized rhs
    nrm = dnrm2(this%ntime, time_function(1), 1)

    if(nrm > CNST(1e-8)) then
      tf_normalized(1:this%ntime) = time_function(1:this%ntime)/nrm
    else
      tf_normalized(1:this%ntime) = time_function(1:this%ntime)
      nrm = M_ONE
    end if

    call bpdn(this%ntime, this%nfreq, this%fourier_matrix, tf_normalized, this%sigma, freq_function, ierr, activesetit = 50)

    SAFE_DEALLOCATE_A(tf_normalized)

    ! scale by the missing factors
    freq_function(1:this%nfreq) = nrm/(this%dfreq*M_TWO/M_PI)*freq_function(1:this%nfreq)

    if(ierr < 0) then
      message(1) = 'The Basis Pursuit Denoising process failed to converge.'
      call messages_warning(1)
    end if

    POP_SUB(compressed_sensing_spectral_analysis)    
  end subroutine compressed_sensing_spectral_analysis

end module compressed_sensing_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
