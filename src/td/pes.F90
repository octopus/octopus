!! Copyright (C) 2006-2011 M. Marques, U. De Giovannini
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

module PES_m
  use batch_m
  use cube_function_m
  use cube_m
  use datasets_m
  use density_m
  use derivatives_m
  use fft_m
  use geometry_m
  use global_m
  use grid_m
  use hamiltonian_m
  use index_m
  use io_binary_m
  use io_function_m
  use io_m
  use lasers_m
  use mesh_m
  use messages_m
  use mpi_m
#if defined(HAVE_NFFT) 
  use nfft_m
#endif
  use output_m
  use parser_m
  use profiling_m
  use qshepmod_m
  use restart_m
  use simul_box_m
  use states_io_m
  use states_m
  use string_m
  use system_m
  use tdpsf_m
  use unit_m
  use unit_system_m
  use varinfo_m
  use varinfo_m

  implicit none

  integer, parameter ::   &
    FREE           =  1,  &    !< The scattering waves evolve in time as free plane waves
    VOLKOV         =  2,  &    !< The scattering waves evolve with exp(i(p-A(t)/c)^2*dt/2)
    CORRECTED1D    =  3,  &
    EMBEDDING1D    =  4,  &     
    VOLKOV_CORRECTED= 5

  integer, parameter ::       &
    PW_MAP_INTEGRAL    =  1,  &    !< projection on outgoing waves by direct integration
    PW_MAP_FFT         =  2,  &    !< FFT on outgoing waves (1D only)
    PW_MAP_BARE_FFT    =  3,  &    !< FFT - normally from fftw3
    PW_MAP_TDPSF       =  4,  &    !< time-dependent phase-space filter
    PW_MAP_NFFT        =  5        !< non-equispaced fft (NFFT)

  integer, parameter ::      &
    M_SIN2            =  1,  &  
    M_STEP            =  2,  & 
    M_ERF             =  3    

  integer, parameter ::       &
    MODE_MASK         =   1,  &  
    MODE_BACKACTION   =   2,  &  
    MODE_PASSIVE      =   3,  &
    MODE_PSF          =   4

  integer, parameter ::      &
    IN                =  1,  &  
    OUT               =  2

  type PES_rc_t
    integer          :: npoints   					!< how many points we store the wf
    integer, pointer :: points(:) 					!< which points to use
    character(len=30), pointer :: filenames(:)  !< filenames
    CMPLX, pointer :: wf(:,:,:,:,:)
    integer, pointer :: rankmin(:)  				   !<partition of the mesh containing the points
  end type PES_rc_t

  type PES_mask_t


    CMPLX, pointer :: k(:,:,:,:,:,:) => NULL() !< The states in momentum space

    ! Some mesh-related stuff
    integer          :: ll(MAX_DIM)            !< the size of the square mesh
    integer          :: np                     !< number of mesh points associated with the mesh
                                               !< (either mesh%np or mesh%np_global)
    FLOAT            :: spacing(MAX_DIM)       !< the spacing
    integer, pointer :: Lxyz_inv(:,:,:)  => NULL()    !< return a point on the main mesh from xyz on the mask square mesh
    type(mesh_t), pointer  :: mesh             !< a pointer to the mesh
    type(cube_t)     :: cube                   !< the cubic mesh

    FLOAT, pointer :: ext_pot(:,:) => NULL()   !< external time-dependent potential i.e. the lasers

    FLOAT, pointer :: M(:,:,:)  => NULL()      !< the mask on a cubic mesh containing the simulation box
    type(cube_function_t) :: cM                !< the mask cube function
    FLOAT, pointer :: mask_R(:)  => NULL()     !< the mask inner (component 1) and outer (component 2) radius
    integer        :: shape                    !< which mask function?

    FLOAT, pointer :: Lk(:) => NULL()          !< associate a k value to an cube index
                                               !< we implicitly assume k to be the same for all directions

    integer          :: resample_lev           !< resampling level
    integer          :: enlarge                !< Fourier space enlargement
    integer          :: enlarge_nfft           !< NFFT space enlargement
    integer          :: llr(MAX_DIM)           !< the size of the rescaled cubic mesh
       
    FLOAT :: energyMax 
    FLOAT :: energyStep 

    integer :: sw_evolve             !< choose the time propagator for the continuum wfs
    logical :: back_action           !< whether to enable back action from B to A
    logical :: add_psia              !< add the contribution of Psi_A in the buffer region to the output
    logical :: interpolate_out       !< whether to apply interpolation on the output files

    integer :: mode                  !< calculation mode
    integer :: pw_map_how            !< how to perform projection on plane waves

    type(fft_t)    :: fft            !< FFT plan
#if defined(HAVE_NFFT) 
    type(nfft_t)   :: nfft           !< NFFT plan
#endif

    type(tdpsf_t) :: psf             !< Phase-space filter struct reference


  end type PES_mask_t




  type PES_t
    logical :: calc_rc
    type(PES_rc_t) :: rc

    logical :: calc_mask
    type(PES_mask_t) :: mask
    
  end type PES_t

  integer, parameter ::     &
    PHOTOELECTRON_NONE = 0, &
    PHOTOELECTRON_RC   = 2, &
    PHOTOELECTRON_MASK = 4

contains

  ! ---------------------------------------------------------
  subroutine PES_init(pes, mesh, sb, st, save_iter,hm, max_iter,dt,sys)
    type(pes_t),         intent(out)   :: pes
    type(mesh_t),        intent(inout) :: mesh
    type(simul_box_t),   intent(in)    :: sb
    type(states_t),      intent(in)    :: st
    integer,             intent(in)    :: save_iter
    type(hamiltonian_t), intent(in)    :: hm
    integer,             intent(in)    :: max_iter
    FLOAT,               intent(in)    :: dt
    type(system_t),      intent(in)    :: sys

    character(len=50)    :: str
    integer :: photoelectron_flags

    PUSH_SUB(PES_init)

    call messages_obsolete_variable('CalcPES_rc', 'PhotoElectronSpectrum')
    call messages_obsolete_variable('CalcPES_mask', 'PhotoElectronSpectrum')

    !%Variable PhotoElectronSpectrum
    !%Type flag
    !%Default no
    !%Section Time-Dependent::PhotoElectronSpectrum
    !%Description
    !%This variable controls the method used for the calculation of
    !%the photoelectron spectrum. You can specify more than one value
    !%by giving them as a sum, for example:
    !% <tt>PhotoElectronSpectrum = pes_rc + pes_mask</tt>
    !%Option none 0
    !% The photoelectron spectrum is not calculated. This is the default.
    !%Option pes_rc 2
    !% Store the wavefunctions at specific points in order to 
    !% calculate the photoelectron spectrum at a point far in the box as proposed in 
    !% A. Pohl, P.-G. Reinhard, and E. Suraud, <i>Phys. Rev. Lett.</i> <b>84</b>, 5090 (2000).
    !%Option pes_mask 4
    !% Calculate the photo-electron spectrum using the mask method.
    !% (D. Varsano, PhD thesis, page 159 (2006) http://nano-bio.ehu.es/files/varsano_phd.pdf).
    !%End

    call parse_integer(datasets_check('PhotoElectronSpectrum'), PHOTOELECTRON_NONE, photoelectron_flags)
    if(.not.varinfo_valid_option('PhotoElectronSpectrum', photoelectron_flags, is_flag = .true.)) then
      call input_error('PhotoElectronSpectrum')
    end if
    
    pes%calc_rc = iand(photoelectron_flags, PHOTOELECTRON_RC) /= 0
    pes%calc_mask = iand(photoelectron_flags, PHOTOELECTRON_MASK) /= 0

    !Header Photoelectron info
    if(pes%calc_rc .or. pes%calc_mask) then 
      write(str, '(a,i5)') 'Photoelectron'
      call messages_print_stress(stdout, trim(str))
    end if 

    
    if(pes%calc_rc) call PES_rc_init(pes%rc, mesh, st, save_iter)
    if(pes%calc_mask) call PES_mask_init(pes%mask, mesh, sb, st,hm,max_iter,dt)


    !Footer Photoelectron info
    if(pes%calc_rc .or. pes%calc_mask) then 
      call messages_print_stress(stdout)
    end if 

    POP_SUB(PES_init)
  end subroutine PES_init


  ! ---------------------------------------------------------
  subroutine PES_end(pes)
    type(PES_t), intent(inout) :: pes

    PUSH_SUB(PES_end)

    if(pes%calc_rc)   call PES_rc_end  (pes%rc)
    if(pes%calc_mask) call PES_mask_end(pes%mask)

    POP_SUB(PES_end)
  end subroutine PES_end


  ! ---------------------------------------------------------
  subroutine PES_calc(pes, mesh, st, ii, dt, mask,hm,geo,iter)
    type(PES_t),         intent(inout) :: pes
    type(mesh_t),        intent(in)    :: mesh
    type(states_t),      intent(inout) :: st
    FLOAT,               intent(in)    :: dt
    FLOAT,               intent(in)    :: mask(:)
    integer,             intent(in)    :: ii
    integer,             intent(in)    :: iter
    type(hamiltonian_t), intent(in)    :: hm
    type(geometry_t),    intent(in)    :: geo

    PUSH_SUB(PES_calc)

    if(pes%calc_rc)   call PES_rc_calc  (pes%rc, st, mesh, ii)
    if(pes%calc_mask) call PES_mask_calc(pes%mask, mesh, st, dt, mask, hm, geo, iter)

    POP_SUB(PES_calc)
  end subroutine PES_calc


  ! ---------------------------------------------------------
  subroutine PES_output(pes, mesh, st, iter, outp, dt, gr, geo)
    type(PES_t),      intent(inout) :: pes
    type(mesh_t),     intent(in)    :: mesh
    type(states_t),   intent(in)    :: st
    integer,          intent(in)    :: iter
    type(output_t),   intent(in)    :: outp
    FLOAT,            intent(in)    :: dt
    type(grid_t),     intent(inout) :: gr
    type(geometry_t), intent(in)    :: geo



    PUSH_SUB(PES_output)
    
    if(mpi_grp_is_root(mpi_world)) then
      if(pes%calc_rc)   call PES_rc_output   (pes%rc, st, iter,outp%iter, dt)
    endif

    if(pes%calc_mask) call PES_mask_output (pes%mask, mesh, st,outp, "td.general/PESM",gr, geo,iter)


    POP_SUB(PES_output)
  end subroutine PES_output

  ! ---------------------------------------------------------
  subroutine PES_restart_write(pes, mesh, st)
    type(PES_t),    intent(in) :: pes
    type(mesh_t),   intent(in) :: mesh
    type(states_t), intent(in) :: st

    PUSH_SUB(PES_restart_write)

      if(pes%calc_mask) call PES_mask_restart_write (pes%mask, mesh, st)

    POP_SUB(PES_restart_write)
  end subroutine PES_restart_write

  ! ---------------------------------------------------------
  subroutine PES_restart_read(pes, mesh, st)
    type(PES_t),    intent(inout) :: pes
    type(mesh_t),   intent(in)    :: mesh
    type(states_t), intent(inout) :: st

    PUSH_SUB(PES_restart_read)

      if(pes%calc_mask) call PES_mask_restart_read (pes%mask, mesh, st)


    POP_SUB(PES_restart_read)
  end subroutine PES_restart_read


  ! ---------------------------------------------------------
  subroutine PES_init_write(pes, mesh, st)
    type(PES_t),    intent(in)  :: pes
    type(mesh_t),   intent(in)  :: mesh
    type(states_t), intent(in)  :: st


    PUSH_SUB(PES_init_write)

    if(mpi_grp_is_root(mpi_world)) then

      if(pes%calc_rc)   call PES_rc_init_write (pes%rc, mesh, st)

    endif

    POP_SUB(PES_init_write)
  end subroutine PES_init_write


#include "pes_rc_inc.F90"
#include "pes_mask_inc.F90"
#include "pes_mask_out_inc.F90"

end module PES_m


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
