!! Copyright (C) 2011 U. De Giovannini
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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!

#include "global.h"

program photoelectron_spectrum
  use command_line_oct_m
  use global_oct_m
  use kpoints_oct_m
  use io_function_oct_m
  use io_oct_m
  use ions_oct_m
  use messages_oct_m
  use multicomm_oct_m
  use namespace_oct_m
  use parser_oct_m
  use pes_mask_oct_m
  use pes_flux_oct_m
  use pes_out_oct_m
  use profiling_oct_m
  use restart_oct_m
  use simul_box_oct_m
  use sort_oct_m
  use space_oct_m
  use string_oct_m
  use states_elec_oct_m
  use states_elec_dim_oct_m
  use symmetries_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use utils_oct_m
  use mpi_oct_m
  
  implicit none

  type pesoutput_t
    logical             :: what(MAX_OUTPUT_TYPES)
    integer(8)          :: how(0:MAX_OUTPUT_TYPES)
    integer             :: output_interval(0:MAX_OUTPUT_TYPES)   
    FLOAT               :: pol(3)
    FLOAT               :: pvec(3)        
  end type pesoutput_t  

  integer              :: ierr, integrate
  integer              :: dim, dir, idim, pdim
  logical              :: what(MAX_OUTPUT_TYPES)
  integer              :: llp(3), llpp(3)  !< The size of the p-point cubic grids 
  FLOAT                :: Emax, Emin, Estep, uEstep,uEspan(2), pol(3)
  FLOAT                :: uThstep, uThspan(2), uPhstep, uPhspan(2), pvec(3)
  FLOAT                :: center(3)
  FLOAT, allocatable   :: Lg(:,:), RR(:)
  FLOAT, allocatable   :: pmesh(:,:,:,:)   !< The final momentum-space (p) mesh 
  integer, allocatable :: Lp(:,:,:,:,:)    !< An index mapping from g- and k-point mesh to p-mesh
  logical              :: need_pmesh, resolve_states
  integer              :: ii, i1,i2,i3, idxZero(1:3), st_range(2)
  type(block_t)        :: blk  
  
  type(space_t)        :: space
  type(ions_t),     pointer :: ions
  type(simul_box_t)    :: sb
  type(states_elec_t)  :: st
  type(symmetries_t)   :: symm
  type(kpoints_t)      :: kpoints
  type(restart_t)      :: restart
  
  character(len=512)   :: filename

  logical              :: have_zweight_path, use_zweight_path 
  integer              :: krng(2), nkpt, kpth_dir !< Kpoint range for zero-weight path

  type(pesoutput_t)    :: pesout
    
  integer              :: ist, ispin  
  FLOAT, pointer       :: pesP_out(:,:,:)    
  FLOAT, allocatable, target :: pesP(:,:,:,:)     !< The momentum-resolved photoelectron spectrum
  FLOAT, allocatable   :: Ekin(:,:,:)    
  
  type(pes_flux_t)     :: pflux
  integer              :: pes_method, option 

  type(multicomm_t)    :: mc
  integer              :: index_range(4)

  call getopt_init(ierr)
  if(ierr /= 0) then
    message(1) = "Your Fortran compiler doesn't support command-line arguments;"
    message(2) = "the oct-photoelectron-spectrum command is not available."
    call messages_fatal(2)
  end if


  call global_init(is_serial = .true.)

  call parser_init()

  call messages_init()  
  call io_init()

  !* In order to initialize k-points
  call unit_system_init(global_namespace)
  
  call space_init(space, global_namespace)
  ions => ions_t(global_namespace)
  call simul_box_init(sb, global_namespace, ions, space)
  call symmetries_init(symm, global_namespace, ions, space)

  ! we need k-points for periodic systems
  call kpoints_init(kpoints, global_namespace, symm, space%dim, space%periodic_dim, ions%latt)
  call states_elec_init(st, global_namespace, space, ions%val_charge(), kpoints)
  !*

  !Initialize variables
  llp(:) = 1 
  llpp(:) = 1

  need_pmesh = .false. 
  dim    = space%dim   ! The dimensionality dim = [1,2,3]
  pdim   = space%periodic_dim

  call messages_print_stress(stdout,"Postprocessing")  
  
  !Figure out which method has been used to calculate the photoelectron data  
  call parse_variable(global_namespace, 'PhotoElectronSpectrum', OPTION__PHOTOELECTRONSPECTRUM__NONE, pes_method)
  
  select case (pes_method)
  case (OPTION__PHOTOELECTRONSPECTRUM__PES_MASK)
    call messages_write('Will process mask-method data.')
    call messages_new_line()  
    
    ! Note that Lg(:,:) is allocated inside pes_mask_read_info
    call pes_mask_read_info("td.general/", global_namespace, dim, Emax, Estep, llp(:), Lg, RR)
    ! Keep a copy the original dimensions vector
    llpp(1:dim) = llp(1:dim) 

    call messages_write('Read PES_MASK info file.')
    call messages_info()
    
    need_pmesh = space%is_periodic() 
    
    
  case (OPTION__PHOTOELECTRONSPECTRUM__PES_FLUX)
    call messages_write('Will process flux-method data.')
    call messages_new_line()
    call messages_info()
    
    option = OPTION__PES_FLUX_SHAPE__SPH
    if(dim <= 2) option = OPTION__PES_FLUX_SHAPE__CUB
    if (space%is_periodic()) option = OPTION__PES_FLUX_SHAPE__PLN
    
    call parse_variable(global_namespace, 'PES_Flux_Shape', option, pflux%surf_shape)
    call pes_flux_reciprocal_mesh_gen(pflux, global_namespace, space, sb, st, kpoints, 0, post = .true.)
    
    llpp(1:dim) = pflux%ll(1:dim)
    need_pmesh = .true.
    
  
  case (OPTION__PHOTOELECTRONSPECTRUM__PES_SPM)
    call messages_not_implemented('Postprocessing SPM data.')  
    call messages_fatal()

  case default 
    call messages_write('Could not find any photoelectron data')
    call messages_fatal()
      
  end select


  !set default values


  integrate = INTEGRATE_NONE
  uEstep  = -1
  uEspan  = (/-1,-1/)
  uThstep = -1
  uThspan = (/-1,-1/)
  uPhstep = -1
  uPhspan = (/-1,-1/)
  center  = (/0,0,0/)
  pvec    = (/1,0,0/)

  what(OPTION__PHOTOELECTRONSPECTRUMOUTPUT__ENERGY_TOT) = .true.
  pesout%pvec = pvec

  have_zweight_path = kpoints%have_zero_weight_path()
  use_zweight_path  = have_zweight_path

  call get_laser_polarization(pol)
  ! if there is no laser set pol along the z-axis
  if (sum(pol(1:3)**2) <= M_EPSILON) pol = (/0,0,1/) 

  ! more defaults
  if(space%is_periodic()) then
    
    if (dim == 2) then
      ! write the velocity map on plane pz=0 as it contains all the informations
      pol  = (/0,1,0/) 
      pvec = (/0,0,1/)
      
      what(OPTION__PHOTOELECTRONSPECTRUMOUTPUT__VELOCITY_MAP_CUT) = .true.
      pesout%pol  = (/0,1,0/) 
      pesout%pvec = (/0,0,1/)
    end if
    
    if (dim == 3) then 
      ! write the full ARPES in vtk format (this could be a big file)
      pol = (/0,0,1/) 
      
      what(OPTION__PHOTOELECTRONSPECTRUMOUTPUT__ARPES) = .true.
      pesout%pol = (/0,0,1/) 
    end if
    
    if (have_zweight_path) then
      ! In this case the output is well defined only on a path in reciprocal space
      ! so we are going to have only a 2D slice regardless of dim=2 or 3 
      pol  = (/0,0,1/) 
      pvec = (/0,1,0/)

      what(OPTION__PHOTOELECTRONSPECTRUMOUTPUT__ARPES_CUT) = .true.
      pesout%pol  = (/0,0,1/) 
      pesout%pvec = (/0,1,0/)
    end if 
    
  end if

  
  ! Intialize mc 
  call multicomm_init(mc, global_namespace, mpi_world, P_STRATEGY_KPOINTS, P_STRATEGY_KPOINTS, &
                      mpi_world%size, index_range, (/ 5000, 1, 1, 1 /))
  index_range(:) = 100000        
  
  
  call restart_module_init(global_namespace)
  call restart_init(restart, global_namespace, RESTART_TD, RESTART_TYPE_LOAD, mc, ierr)
  if(ierr /= 0) then
    message(1) = "Unable to read time-dependent restart information."
    call messages_fatal(1)
  end if
  
  !%Variable PhotoelectronSpectrumOutput
  !%Type block
  !%Default none
  !%Section Utilities::oct-photoelectron_spectrum
  !%Description
  !% Specifies what to output extracting the photoelectron cross-section informations. 
  !% When we use polar coordinates the zenith axis is set by vec (default is the first 
  !% laser field polarization vector), theta is the inclination angle measured from 
  !% vec (from 0 to \pi), and phi is the azimuthal angle on a plane perpendicular to 
  !% vec (from 0 to 2\pi).
  !% Each option must be in a separate row. Optionally individual output formats can be defined
  !% for each row or they can be read separately from <tt>OutputFormat</tt> variable
  !% in the input file.
  !%
  !% Example (minimal):
  !% <br><br><tt>%PhotoelectronSpectrumOutput
  !% <br>&nbsp;&nbsp;energy_tot
  !% <br>&nbsp;&nbsp;velocity_map
  !% <br>%<br></tt>
  !%
  !% Example (with OutputFormat):
  !% <br><br><tt>%PhotoelectronSpectrumOutput
  !% <br>&nbsp;&nbsp;arpes        | vtk
  !% <br>&nbsp;&nbsp;velocity_map | ncdf
  !% <br>%<br></tt>
  !%
  !%Option energy_tot 1
  !% Output the energy-resolved photoelectron spectrum: E.
  !%Option energy_angle 2
  !% Output the energy and angle resolved spectrum: (theta, E)
  !% The result is integrated over phi.
  !%Option velocity_map_cut 3
  !% Velocity map on a plane orthogonal to pvec: (px, py). The allowed cutting planes 
  !% (pvec) can only be parallel to the x,y,z=0 planes. 
  !% Space is oriented so that the z-axis is along vec. Supports the -I option.
  !%Option energy_xy 4
  !% Angle and energy-resolved spectrum on the inclination plane: (Ex, Ey).
  !% The result is integrated over ph;
  !%Option energy_th_ph 5
  !% Ionization probability integrated on spherical cuts: (theta, phi).
  !%Option velocity_map 6
  !% Full momentum-resolved ionization probability: (px, py, pz).     
  !% The output format can be controlled with <tt>OutputHow</tt> and can be vtk, ncdf or ascii.  
  !%Option arpes 7
  !% Full ARPES for semi-periodic systems (vtk).
  !%Option arpes_cut 8
  !% ARPES cut on a plane following a zero-weight path in reciprocal space.
  !%End
  pesout%what = what
  call io_function_read_what_how_when(global_namespace, space, pesout%what, pesout%how, pesout%output_interval, &
    what_tag_in = 'PhotoelectronSpectrumOutput',  ignore_error = .true.)
  
  ! TODO: I think it would be better to move these options in the
  ! input file to have more flexibility to combine and to keep
  ! track of them. UDG
  ! Read options from command line
  call getopt_photoelectron_spectrum(uEstep, uEspan,&
                                     uThstep, uThspan, uPhstep, &
                                     uPhspan, pol, center, pvec, integrate)
                                     
  call getopt_end()
                                       

  !! set user values
  if(uEstep >  0 .and. uEstep > Estep)    Estep = uEstep
  if(uEspan(1) > 0 ) Emin = uEspan(1)
  if(uEspan(2) > 0 ) Emax = uEspan(2)
  
  
  
  !%Variable PhotoelectronSpectrumResolveStates
  !%Type block
  !%Default 
  !%Section Utilities::oct-photoelectron_spectrum
  !%Description
  !% If </tt>yes</tt> calculate the photoelectron spectrum resolved in each K.S. state. 
  !% Optionally a range of states can be given as two slot block where the 
  !% first slot is the lower state index and the second is the highest one. 
  !% For example to calculate the spectra from state i to state j:
  !%
  !% <tt>%PhotoelectronSpectrumResolveStates
  !% <br> i | j
  !% <br>%</tt>
  !%End
  st_range(1:2) = (/1, st%nst/)
  resolve_states = .false.
  if(parse_block(global_namespace, 'PhotoelectronSpectrumResolveStates', blk) == 0) then
    if(parse_block_cols(blk,0) < 2) call messages_input_error(global_namespace, 'PhotoelectronSpectrumResolveStates')
    do idim = 1, 2
      call parse_block_integer(blk, 0, idim - 1, st_range(idim))
    end do
    call parse_block_end(blk)
    if (abs(st_range(2)-st_range(1)) > 0)resolve_states = .true.    
  else
    call parse_variable(global_namespace, 'PhotoelectronSpectrumResolveStates', .false., resolve_states)
  end if
  
  
  krng(1) = 1
  krng(2) = kpoints_number(kpoints)
  
  if (have_zweight_path) then 
    
    if (pesout%what(OPTION__PHOTOELECTRONSPECTRUMOUTPUT__ARPES_CUT)) then
      !Use the path only when asked for ARPES on a cutting curve in reciprocal space(the path)
      use_zweight_path  = .true.
    
      krng(1) = kpoints_number(kpoints) - kpoints%nik_skip  + 1
    
      call messages_print_stress(stdout, "Kpoint selection")
      write(message(1), '(a)') 'Will use a zero-weight path in reciprocal space with the following points'
      call messages_info(1)

      kpth_dir = 1
      pvec     = (/0,1,0/)
    

    else
      use_zweight_path  = .false.
      krng(2) = kpoints_number(kpoints) - kpoints%nik_skip

    end if
    
  end if

  call write_kpoints_info(kpoints, krng(1), krng(2))    
  call messages_print_stress(stdout)

  
  nkpt = krng(2) - krng(1) + 1

  
 
  if (use_zweight_path) then
    llp(1:dim) = llpp(1:dim)
    llp(kpth_dir) = llpp(kpth_dir) * nkpt    
  else
    llp(1:dim) = llpp(1:dim)  
  endif  
  
  if (debug%info) then
    write(message(1),'(a,i4,i4,i4)') 'Debug :  llp = ', llp(1:3) 
    write(message(2),'(a,i4,i4,i4)') 'Debug : llpp = ', llpp(1:3) 
    call messages_info(2)
  end if
  
  SAFE_ALLOCATE(pmesh(1:llp(1), 1:llp(2), 1:llp(3), 1:3 + 1))
  SAFE_ALLOCATE( pesP(1:llp(1), 1:llp(2), 1:llp(3), 1:st%d%nspin))

  select case (pes_method)
  case (OPTION__PHOTOELECTRONSPECTRUM__PES_MASK)
    SAFE_ALLOCATE(Lp(1:llpp(1), 1:llpp(2), 1:llpp(3), krng(1):krng(2), 1:3))
    call pes_mask_pmesh(global_namespace, dim, kpoints, llpp, Lg, pmesh, idxZero, krng, Lp)

  case (OPTION__PHOTOELECTRONSPECTRUM__PES_FLUX)
    ! Lp is allocated inside pes_flux_pmesh to comply with the 
    ! declinations of the different surfaces
    SAFE_ALLOCATE(Ekin(1:llp(1), 1:llp(2), 1:llp(3)))
    Ekin = M_ZERO
    call pes_flux_pmesh(pflux, global_namespace, dim, kpoints, llpp, pmesh, idxZero, krng, Lp, Ekin)
  end select
   
  
  if (.not. need_pmesh) then
    ! There is no need to use pmesh we just need to sort Lg in place
    ! in order to have a coordinate ordering coherent with pesP
    ! NOTE: this works only for the mask_method since Lg is well-defined  
    ! only in this case
    do idim = 1, dim
      call sort(Lg(1:llp(idim), idim)) 
    end do  
  end if  


  call messages_write('Read PES restart files.')
  call messages_info()



  call unit_system_init(global_namespace)
 
  write(message(1),'(a,f10.2,a2,f10.2,a2,f10.2,a1)') &
                   "Zenith axis: (",pol(1),", ",pol(2),", ",pol(3),")"
  call messages_info(1)


  ! Convert the grid units
  if (need_pmesh) then    
    do ii = 1,3
      do i3 = 1, llp(3)
        do i2 = 1, llp(2)
          do i1 = 1, llp(1)
            pmesh(i1, i2, i3, ii) = units_from_atomic(sqrt(units_out%energy), pmesh(i1, i2, i3, ii))
          end do
        end do
      end do
    end do
  end if

  if (resolve_states) then
    do ist = st_range(1), st_range(2)

      select case (pes_method)
      case (OPTION__PHOTOELECTRONSPECTRUM__PES_MASK)
        call pes_mask_map_from_states(restart, st, llpp, pesP, krng, Lp, ist)
      case (OPTION__PHOTOELECTRONSPECTRUM__PES_FLUX)
        call pes_flux_map_from_states(pflux, restart, st, llpp, pesP, krng, Lp, ist)      
      end select

      call output_spin_pes()
    end do

  else
    ! Read the data
    ist = 0 

    select case (pes_method)
    case (OPTION__PHOTOELECTRONSPECTRUM__PES_MASK)
      call pes_mask_map_from_states(restart, st, llpp, pesP, krng, Lp)
    case (OPTION__PHOTOELECTRONSPECTRUM__PES_FLUX)
      call pes_flux_map_from_states(pflux, restart, st, llpp, pesP, krng, Lp)      
    end select

    call output_spin_pes()

  end if



  write(message(1), '(a)') 'Done'
  call messages_info(1)

  call messages_print_stress(stdout)

  call restart_end(restart)    

  call states_elec_end(st)

  SAFE_DEALLOCATE_P(ions)
  call kpoints_end(kpoints)
  call simul_box_end(sb)

  call io_end()
  call messages_end()

  call parser_end()
  call global_end()
  
  SAFE_DEALLOCATE_A(pesP)    
  SAFE_DEALLOCATE_A(pmesh)
  SAFE_DEALLOCATE_A(Lp)
  SAFE_DEALLOCATE_A(Ekin)
  if (.not. need_pmesh .or. pes_method == OPTION__PHOTOELECTRONSPECTRUM__PES_MASK) then
    SAFE_DEALLOCATE_A(Lg)
  end if
  contains
    
    function outfile(name, ist, ispin, extension) result(fname)    
      character(len=*),   intent(in)   :: name
      integer,            intent(in)   :: ist
      integer,            intent(in)   :: ispin
      character(len=*), optional, intent(in)    :: extension
      character(len=512)               :: fname
            
      if (ist == 0) then 
        fname = trim(name) 
      else
        write(fname, '(a,a,i4.4)') trim(name), '-st', ist
      end if 

      if (ispin >0)  then
        write(fname, '(a,a,i1)') trim(fname),'-sp', ispin
      end if 
      
      if (present(extension) .and. len(extension)>0) then
        fname = trim(fname)//'.'//trim(extension)
      end if
          
    end function outfile  
          
    subroutine output_spin_pes()

      ! Write total quantities (summed over spin) 
      ispin = 0
      if (st%d%ispin == UNPOLARIZED) then
        pesP_out => pesP(:,:,:,1)
        call output_pes()
      else 
        ! Write total quantities (summed over spin) 
        SAFE_ALLOCATE(pesP_out(1:llp(1), 1:llp(2), 1:llp(3)))
        pesP_out(:,:,:) = pesP(:,:,:,1) + pesP(:,:,:,2)
    
        call output_pes()
        
        SAFE_DEALLOCATE_P(pesP_out)      
    
        ! spin-resolved 
        do ispin = 1, st%d%nspin
          pesP_out => pesP(:,:,:,ispin)
          call output_pes()    
        end do
      end if      
    end subroutine output_spin_pes
    
    subroutine output_pes()
      
      ! choose what to calculate
      ! these functions are defined in pes_mask_out_inc.F90

      if (st%d%ispin /= UNPOLARIZED .or. ist>0) call messages_print_stress(stdout)
      
      if (ist > 0 ) then 
        write(message(1), '(a,i4)') 'State = ', ist
        call messages_info(1)
      end if
      
      if (st%d%ispin /= UNPOLARIZED) then
        if (ispin > 0 ) then 
          write(message(1), '(a,i1)') 'Spin component= ', ispin
          call messages_info(1)
        else 
          write(message(1), '(a)') 'Spinless'
          call messages_info(1)
        end if
      end if
      
      if (pesout%what(OPTION__PHOTOELECTRONSPECTRUMOUTPUT__ENERGY_TOT)) then
        call messages_print_stress(stdout, "Energy-resolved PES")

        select case (pes_method)
        case (OPTION__PHOTOELECTRONSPECTRUM__PES_MASK)
          call pes_mask_output_power_totalM(pesP_out,outfile('./PES_energy',ist, ispin, 'sum'), &
                                            global_namespace, Lg, llp, dim, Emax, Estep, interpolate = .true.)
        case (OPTION__PHOTOELECTRONSPECTRUM__PES_FLUX)
          call pes_flux_out_energy(pflux, pesP_out, outfile('./PES_energy',ist, ispin, 'sum'), global_namespace, llp, Ekin, dim)
        end select 
        
      end if
      
      if (pesout%what(OPTION__PHOTOELECTRONSPECTRUMOUTPUT__ENERGY_ANGLE)) then
        call messages_print_stress(stdout, "Angle- and energy-resolved PES")
        
        select case (pes_method)
        case (OPTION__PHOTOELECTRONSPECTRUM__PES_MASK)
          call pes_mask_output_ar_polar_M(pesP_out,outfile('./PES_angle_energy',ist, ispin, 'map'), &
                                          global_namespace, Lg, llp, dim, pol, Emax, Estep)
        case (OPTION__PHOTOELECTRONSPECTRUM__PES_FLUX)                                     
          call messages_not_implemented("Angle- and energy-resolved PES for the flux method") 
        end select                       
                                          
                                          
      end if

      if (pesout%what(OPTION__PHOTOELECTRONSPECTRUMOUTPUT__VELOCITY_MAP_CUT)) then
        call messages_print_stress(stdout, "Velocity map on a plane")
        dir = -1
        if(sum((pvec-(/1 ,0 ,0/))**2)  <= M_EPSILON  )  dir = 1
        if(sum((pvec-(/0 ,1 ,0/))**2)  <= M_EPSILON  )  dir = 2
        if(sum((pvec-(/0 ,0 ,1/))**2)  <= M_EPSILON  )  dir = 3

        if (use_zweight_path) then
          filename = outfile('PES_velocity_map', ist, ispin, 'path')
        else
          filename = outfile('PES_velocity_map', ist, ispin, 'p'//index2axis(dir)//'=0')
        end if

        if (dir == -1) then
            write(message(1), '(a)') 'Unrecognized plane. Use -u to change.'
            call messages_fatal(1)
          else
            write(message(1), '(a)') 'Save velocity map on plane: '//index2axis(dir)//" = 0"
            call messages_info(1)
        end if

        if(integrate /= INTEGRATE_NONE) then
          write(message(1), '(a)') 'Integrate on: '//index2var(integrate)
          call messages_info(1)
          filename = trim(filename)//'.i_'//trim(index2var(integrate))
        end if

        if (need_pmesh) then
          call pes_out_velocity_map_cut(global_namespace, pesP_out, filename, llp, dim, pol, dir, integrate, &
                                             pos = idxZero, pmesh = pmesh)
        else
          call pes_out_velocity_map_cut(global_namespace, pesP_out, filename, llp, dim, pol, dir, integrate, &
                                             pos = idxZero, Lk = Lg)
        end if
      end if

      if (pesout%what(OPTION__PHOTOELECTRONSPECTRUMOUTPUT__ENERGY_XY)) then
        call messages_print_stress(stdout, "Angle and energy-resolved on a plane")
        if(uEstep >  0 .and. uEstep > Estep) then
          Estep = uEstep
        else
          Estep = Emax/size(Lg,1)
        end if

        select case (pes_method)
        case (OPTION__PHOTOELECTRONSPECTRUM__PES_MASK)
          call pes_mask_output_ar_plane_M(pesP_out,outfile('./PES_energy', ist, ispin, 'map'), &
                                          global_namespace, Lg, llp, dim, pol, Emax, Estep)
        case (OPTION__PHOTOELECTRONSPECTRUM__PES_FLUX)                                     
          call messages_not_implemented("Angle and energy-resolved on a plane for the flux method") 
        end select   
                                        
      end if

      if (pesout%what(OPTION__PHOTOELECTRONSPECTRUMOUTPUT__ENERGY_TH_PH)) then
        call messages_print_stress(stdout, "PES on spherical cuts")

        write(message(1), '(a,es19.12,a2,es19.12,2x,a19)') &
              'Save PES on a spherical cut at E= ',Emin,", ",Emax, &
               str_center('['//trim(units_abbrev(units_out%energy)) // ']', 19)
        call messages_info(1)

        if(uEstep >  0 .and. uEstep > Estep) then
         Estep = uEstep
        else
         Estep = Emax/size(Lg,1)
        end if
        
        select case (pes_method)
        case (OPTION__PHOTOELECTRONSPECTRUM__PES_MASK)
          call pes_mask_output_ar_spherical_cut_M(pesP_out,outfile('./PES_sphere', ist, ispin, 'map'), & 
                                                  global_namespace, Lg, llp, dim, pol, Emin, Emax, Estep)

        case (OPTION__PHOTOELECTRONSPECTRUM__PES_FLUX)                                     
          call messages_not_implemented("PES on spherical cuts for the flux method") 
        end select                                          

      end if

      if (pesout%what(OPTION__PHOTOELECTRONSPECTRUMOUTPUT__VELOCITY_MAP)) then
        
        call messages_print_stress(stdout, "Full velocity map")
        
        if ( .not. (bitand(pesout%how(OPTION__PHOTOELECTRONSPECTRUMOUTPUT__VELOCITY_MAP), OPTION__OUTPUTFORMAT__NETCDF) /= 0) .and. &
             .not. (bitand(pesout%how(OPTION__PHOTOELECTRONSPECTRUMOUTPUT__VELOCITY_MAP), OPTION__OUTPUTFORMAT__VTK)    /= 0) .and. &
             .not. (bitand(pesout%how(OPTION__PHOTOELECTRONSPECTRUMOUTPUT__VELOCITY_MAP), OPTION__OUTPUTFORMAT__ASCII)  /= 0) ) then
             message(1) = 'User must specify the format with "OutputFormat".'
             message(2) = 'Available options are: necdf, vtk, ascii.'
             call messages_fatal(2)
             
        end if
        
        filename = outfile('./PES_velocity_map', ist, ispin)
        if (need_pmesh) then
          !force vtk output
!           how = io_function_fill_how("VTK")
          if (bitand(pesout%how(OPTION__PHOTOELECTRONSPECTRUMOUTPUT__VELOCITY_MAP), OPTION__OUTPUTFORMAT__ASCII) /= 0) then
             call pes_flux_out_vmap(pflux, pesP_out, filename, global_namespace, llp, pmesh, space%dim)
          else            
            call pes_out_velocity_map(pesP_out, filename, global_namespace, space, Lg, llp, pesout%how(OPTION__PHOTOELECTRONSPECTRUMOUTPUT__VELOCITY_MAP), pmesh)
          end if
        else
          call pes_out_velocity_map(pesP_out, filename, global_namespace, space, Lg, llp, pesout%how(OPTION__PHOTOELECTRONSPECTRUMOUTPUT__VELOCITY_MAP))
        end if
        
      end if

      if (pesout%what(OPTION__PHOTOELECTRONSPECTRUMOUTPUT__ARPES)) then
        call messages_print_stress(stdout, "ARPES")

        do i3 = 1, llp(3)
          do i2 = 1, llp(2)
            do i1 = 1, llp(1)
              pmesh(i1, i2, i3, dim) = units_from_atomic(units_out%energy, &
                sign(M_ONE,pmesh( i1, i2, i3, dim)) * sum( pmesh(i1, i2, i3, 1:dim)**2 )/M_TWO)
            end do
          end do
        end do

        pesout%how(OPTION__PHOTOELECTRONSPECTRUMOUTPUT__ARPES) = io_function_fill_how("VTK")

        call pes_out_velocity_map(pesP_out, outfile('./PES_ARPES', ist, ispin), &
                                  global_namespace, space, Lg, llp, pesout%how(OPTION__PHOTOELECTRONSPECTRUMOUTPUT__ARPES), pmesh)
      end if
      
      
      if (pesout%what(OPTION__PHOTOELECTRONSPECTRUMOUTPUT__ARPES_CUT)) then
        call messages_print_stress(stdout, "ARPES cut on reciprocal space path")
        
        filename = outfile('./PES_ARPES', ist, ispin, "path")
        call pes_out_arpes_cut(global_namespace, pesP_out, filename, dim, llp, pmesh, Ekin)
        
      end if
      
      
    end subroutine output_pes
    
    
    subroutine write_kpoints_info(kpoints, ikstart, ikend)
      type(kpoints_t),   intent(in) :: kpoints 
      integer,           intent(in) :: ikstart   
      integer,           intent(out):: ikend  

      integer :: ik, idir
      character(len=100) :: str_tmp
            
      PUSH_SUB(write_kpoints_info)

      do ik = ikstart, ikend
        write(message(1),'(i8,1x)') ik
        write(str_tmp,'(f12.6)') kpoints%get_weight(ik)
        message(1) = trim(message(1)) // trim(str_tmp)//' |'
        do idir = 1, kpoints%full%dim
          write(str_tmp,'(f12.6)') kpoints%reduced%red_point(idir, ik)
          message(1) = trim(message(1)) // trim(str_tmp)
        end do
        message(1) = trim(message(1)) //' |'
        call messages_info(1)
      end do
      
      call messages_info(1)
      
      POP_SUB(write_kpoints_info)
    end subroutine write_kpoints_info

    subroutine get_laser_polarization(lPol)
       FLOAT,   intent(out) :: lPol(:) 
       
        type(block_t)       :: blk
        integer             :: no_l
        CMPLX               :: cPol(1:3)
        
        PUSH_SUB(get_laser_polarization)
        
        cPol = M_ZERO

        no_l = 0
        if(parse_block(global_namespace, 'TDExternalFields', blk) == 0) then
          no_l = parse_block_n(blk)

          call parse_block_cmplx(blk, 0, 1, cPol(1))
          call parse_block_cmplx(blk, 0, 2, cPol(2))
          call parse_block_cmplx(blk, 0, 3, cPol(3))


          call parse_block_end(blk)
        end if
        
        lPol(:) = abs(cPol)
        
        if(no_l > 1) then
          message(1)="There is more than one external field. Polarization will be selected"
          message(2)="from the first field. Use -V to change axis."
          call messages_info(2)
        end if

        POP_SUB(get_laser_polarization)
    end subroutine get_laser_polarization

    character(5) pure function index2var(ivar) result(ch)
      integer, intent(in) :: ivar

      select case(ivar)
        case(INTEGRATE_PHI)
          ch = 'phi'
        case(INTEGRATE_THETA)
          ch = 'theta'
        case(INTEGRATE_R)
          ch = 'r'
        case(INTEGRATE_KX)
          ch = 'kx'
        case(INTEGRATE_KY)
          ch = 'ky'
        case(INTEGRATE_KZ)
          ch = 'kz'
        case default
          write(ch,'(i1)') ivar
      end select
    end function index2var
    
    
    



end program photoelectron_spectrum

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
