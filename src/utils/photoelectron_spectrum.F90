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
  use geometry_oct_m
  use global_oct_m
  use grid_oct_m
  use kpoints_oct_m
  use io_function_oct_m
  use io_oct_m
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
  use states_oct_m
  use states_dim_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use utils_oct_m
  
  implicit none

  type pesoutput_t
    integer     :: what  
    FLOAT       :: pol(3)
    FLOAT       :: pvec(3)        
  end type pesoutput_t  

  integer              :: ierr, integrate
  integer              :: dim, dir, idim, pdim, ngpt
  integer(8)           :: how
  integer              :: llp(3), llg(3)   !< The size of the g-point and p-point cubic grids 
  FLOAT                :: Emax, Emin, Estep, uEstep,uEspan(2), pol(3)
  FLOAT                :: uThstep, uThspan(2), uPhstep, uPhspan(2), pvec(3)
  FLOAT                :: center(3)
  FLOAT, pointer       :: Lg(:,:), RR(:)
  FLOAT, allocatable   :: pmesh(:,:,:,:)   !< The final momentum-space (p) mesh 
  integer, pointer     :: Lp(:,:,:,:,:)    !< An index mapping from g- and k-point mesh to p-mesh
  logical              :: need_pmesh, resolve_states
  integer              :: ii, i1,i2,i3, idxZero(1:3), st_range(2)
  type(block_t)        :: blk  
  
  type(space_t)        :: space
  type(geometry_t)     :: geo
  type(simul_box_t)    :: sb
  type(states_t)       :: st
  type(grid_t)         :: gr
  type(restart_t)      :: restart
  
  character(len=512)   :: filename

  logical              :: have_zweight_path 
  integer              :: krng(2), nkpt, kpth_dir !< Kpoint range for zero-weight path

  type(pesoutput_t)    :: pesout
    
  integer              :: ist, ispin  
  FLOAT, pointer       :: pesP_out(:,:,:)    
  FLOAT, allocatable, target :: pesP(:,:,:,:)     !< The momentum-resolved photoelectron spectrum
  FLOAT, allocatable   :: Ekin(:,:,:)    
  
  type(pes_flux_t)     :: pflux
  integer              :: pes_method, option 

  type(multicomm_t)    :: mc
  type(namespace_t) :: namespace
  
  call getopt_init(ierr)
  if(ierr /= 0) then
    message(1) = "Your Fortran compiler doesn't support command-line arguments;"
    message(2) = "the oct-photoelectron-spectrum command is not available."
    call messages_fatal(2)
  end if


  call global_init(is_serial = .true.)

  call parser_init()
  namespace = namespace_t("")

  call messages_init(namespace)  
  call io_init(namespace)

  !* In order to initialize k-points
  call unit_system_init(namespace)
  
  call space_init(space, namespace)
  call geometry_init(geo, namespace, space)
  call simul_box_init(sb, namespace, geo, space)
  gr%sb = sb
  call states_init(st, namespace, gr, geo)
  !*

  !Initialize variables
  llp(:) = 1 
  llg(:) = 1
  ngpt = 1

  need_pmesh = .false. 
  dim    = sb%dim   ! The dimensionality dim = [1,2,3]
  pdim   = sb%periodic_dim

  call messages_print_stress(stdout,"Postprocessing")  
  
  !Figure out wich method has been used to calculate the photoelectron data  
  call parse_variable(namespace, 'PhotoElectronSpectrum', OPTION__PHOTOELECTRONSPECTRUM__NONE, pes_method)
  
  select case (pes_method)
  case (OPTION__PHOTOELECTRONSPECTRUM__PES_MASK)
    call messages_write('Will process mask-method data.')
    call messages_new_line()  
    
    ! Note that Lg(:,:) is allocated inside pes_mask_read_info
    call pes_mask_read_info("td.general/", dim, Emax, Estep, llp(:), Lg, RR)
    ! Keep a copy the original dimensions vector
    ! For periodic systems llg represents the extension on the g-point grid
    llg(1:dim) = llp(1:dim) 

    call messages_write('Read PES_MASK info file.')
    call messages_info()
    
    need_pmesh = simul_box_is_periodic(sb) 
    
    
  case (OPTION__PHOTOELECTRONSPECTRUM__PES_FLUX)
    call messages_write('Will process flux-method data.')
    call messages_new_line()
    call messages_info()
    
    option = OPTION__PES_FLUX_SHAPE__SPH
    if(dim <= 2) option = OPTION__PES_FLUX_SHAPE__CUB
    if (simul_box_is_periodic(sb)) option = OPTION__PES_FLUX_SHAPE__PLN
    
    call parse_variable(namespace, 'PES_Flux_Shape', option, pflux%shape)
    call pes_flux_reciprocal_mesh_gen(pflux, namespace, sb, st, 0, post = .true.)
    
    llg(1:dim) = pflux%ll(1:dim)
    ngpt = pflux%ngpt
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
  uEstep = -1
  uEspan = (/-1,-1/)
  uThstep = -1
  uThspan = (/-1,-1/)
  uPhstep = -1
  uPhspan = (/-1,-1/)
  center = (/0,0,0/)
  pvec = (/1,0,0/)

  pesout%what = OPTION__PHOTOELECTRONSPECTRUMOUTPUT__ENERGY_TOT 
  pesout%pvec = pvec

  have_zweight_path = kpoints_have_zero_weight_path(sb%kpoints)

  call get_laser_polarization(pol)
  ! if there is no laser set pol along the z-axis
  if (sum(pol(1:3)**2) <= M_EPSILON) pol = (/0,0,1/) 

  ! more defaults
  if(simul_box_is_periodic(sb)) then
    
    if (dim == 2) then
      ! write the velocity map on plane pz=0 as it contains all the informations
      pol = (/0,1,0/) 
      pvec = (/0,0,1/)
      
      pesout%what = OPTION__PHOTOELECTRONSPECTRUMOUTPUT__VELOCITY_MAP_CUT 
      pesout%pol = (/0,1,0/) 
      pesout%pvec = (/0,0,1/)
    end if
    
    if (dim == 3) then 
      ! write the full ARPES in vtk format (this could be a big file)
      pol = (/0,0,1/) 
      
      pesout%what = OPTION__PHOTOELECTRONSPECTRUMOUTPUT__ARPES
      pesout%pol = (/0,0,1/) 
    end if
    
    if (have_zweight_path) then
      ! In this case the output is well defined only on a path in reciprocal space
      ! so we are going to have only a 2D slice regardless of dim=2 or 3 
      pol = (/0,0,1/) 
      pvec = (/0,1,0/)

      pesout%what = OPTION__PHOTOELECTRONSPECTRUMOUTPUT__ARPES_CUT
      pesout%pol = (/0,0,1/) 
      pesout%pvec = (/0,1,0/)
    end if 
    
  end if

  
  
  call restart_module_init(namespace)
  call restart_init(restart, namespace, RESTART_TD, RESTART_TYPE_LOAD, mc, ierr)
  if(ierr /= 0) then
    message(1) = "Unable to read time-dependent restart information."
    call messages_fatal(1)
  end if
  
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
  st_range(1:2)=(/1, st%nst/)
  resolve_states = .false.
  if(parse_block(namespace, 'PhotoelectronSpectrumResolveStates', blk) == 0) then
    if(parse_block_cols(blk,0) < 2) call messages_input_error('PhotoelectronSpectrumResolveStates')
    do idim = 1, 2
      call parse_block_integer(blk, 0, idim - 1, st_range(idim))
    end do
    call parse_block_end(blk)
    if (abs(st_range(2)-st_range(1)) > 0)resolve_states = .true.    
  else
    call parse_variable(namespace, 'PhotoelectronSpectrumResolveStates', .false., resolve_states)
  end if
  
  
  krng(1) = 1
  krng(2) =  kpoints_number(sb%kpoints)
  
  if (have_zweight_path) then 
    
    
    krng(1) = kpoints_number(sb%kpoints) - sb%kpoints%nik_skip  + 1
    
    call messages_print_stress(stdout, "Kpoint selection")
    write(message(1), '(a)') 'Will use a zero-weight path in reciprocal space with the following points'
    call messages_info(1)
    ! Figure out the direction of the path - it must be along kx or ky only
!     call get_kpath_direction(sb%kpoints, krng, kpth_dir, pvec)
    kpth_dir = 1
    pvec = (/0,1,0/)
    
    call write_kpoints_info(sb%kpoints, krng(1), krng(2))    
    call messages_print_stress(stdout)
    
  end if
  
  nkpt = krng(2) - krng(1) + 1

  
 
  if (have_zweight_path) then
    llp(1:dim) = llg(1:dim)
    llp(kpth_dir) = llg(kpth_dir) * nkpt    
  else
    llp(1:dim) = llg(1:dim) * sb%kpoints%nik_axis(1:dim)    
  endif  
  llp(1:pdim) = llp(1:pdim) * ngpt
  
  if (debug%info) then
    print *, "llp(:)= ", llp(:) 
    print *, "llg(:)= ", llg(:) 
  end if
  
  SAFE_ALLOCATE(pmesh(1:llp(1),1:llp(2),1:llp(3),1:3 + 1))
  SAFE_ALLOCATE(pesP(1:llp(1),1:llp(2),1:llp(3),1:st%d%nspin))

  select case (pes_method)
  case (OPTION__PHOTOELECTRONSPECTRUM__PES_MASK)
    SAFE_ALLOCATE(Lp(1:llg(1),1:llg(2),1:llg(3),krng(1):krng(2),1:3))
    call pes_mask_pmesh(dim, sb%kpoints, llg, Lg, pmesh, idxZero, krng, Lp)  

  case (OPTION__PHOTOELECTRONSPECTRUM__PES_FLUX)
    ! Lp is allocated inside pes_flux_pmesh to comply with the 
    ! declinations of the different surfaces
    SAFE_ALLOCATE(Ekin(1:llp(1),1:llp(2),1:llp(3)))
    Ekin = M_ZERO
    call pes_flux_pmesh(pflux, dim, sb%kpoints, llg, Lg, pmesh, idxZero, krng, Lp, Ekin)    
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

  !%Variable PhotoelectronSpectrumOutput
  !%Type flag
  !%Default none
  !%Section Utilities::oct-photoelectron_spectrum
  !%Description
  !% Specifies what to output extracting the photoelectron cross-section informations. 
  !% When we use polar coordinates the zenith axis is set by vec (default is the first 
  !% laser field polarization vector), theta is the inclination angle measured from 
  !% vec (from 0 to \pi), and phi is the azimuthal angle on a plane perpendicular to 
  !% vec (from 0 to 2\pi).
  !% Example: <tt>energy_tot + velocity_map</tt>
  !%Option energy_tot bit(1)
  !% Output the energy-resolved photoelectron spectrum: E.
  !%Option energy_angle bit(2)
  !% Output the energy and angle resolved spectrum: (theta, E)
  !% The result is integrated over phi.
  !%Option velocity_map_cut bit(3)
  !% Velocity map on a plane orthogonal to pvec: (px, py). The allowed cutting planes 
  !% (pvec) can only be parallel to the x,y,z=0 planes. 
  !% Space is oriented so that the z-axis is along vec. Supports the -I option.
  !%Option energy_xy bit(4)  
  !% Angle and energy-resolved spectrum on the inclination plane: (Ex, Ey).
  !% The result is integrated over ph;
  !%Option energy_th_ph bit(5)  
  !% Ionization probability integrated on spherical cuts: (theta, phi).
  !%Option velocity_map bit(6)
  !% Full momentum-resolved ionization probability: (px, py, pz).     
  !% The output format can be controlled with OutputHow and can be vtk or ncdf.  
  !%Option arpes bit(7)
  !% Full ARPES for semi-periodic systems (vtk).
  !%Option arpes_cut bit(8)
  !% ARPES cut on a plane following a zero-weight path in reciprocal space.
  !%End
  call parse_variable(namespace, 'PhotoelectronSpectrumOutput', pesout%what, pesout%what)
  
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


  call unit_system_init(namespace)
 
  write(message(1),'(a,f10.2,a2,f10.2,a2,f10.2,a1)') &
                   "Zenith axis: (",pol(1),", ",pol(2),", ",pol(3),")"
  call messages_info(1)


  ! Convert the grid units
  if (need_pmesh) then    
    forall (i1=1:llp(1), i2=1:llp(2), i3=1:llp(3), ii = 1:3)
      pmesh(i1,i2,i3,ii) = units_from_atomic(sqrt(units_out%energy), pmesh(i1,i2,i3,ii))
    end forall
  end if

  if (resolve_states) then
    do ist = st_range(1), st_range(2)
      
      select case (pes_method)
      case (OPTION__PHOTOELECTRONSPECTRUM__PES_MASK)
        call pes_mask_map_from_states(restart, st, llg, pesP, krng, Lp, ist)
      case (OPTION__PHOTOELECTRONSPECTRUM__PES_FLUX)
        call pes_flux_map_from_states(pflux, restart, st, llg, pesP, krng, Lp, ist)      
      end select
        
      call output_spin_pes()
    end do
    
  else
    ! Read the data
    ist = 0 

    select case (pes_method)
    case (OPTION__PHOTOELECTRONSPECTRUM__PES_MASK)
      call pes_mask_map_from_states(restart, st, llg, pesP, krng, Lp)
    case (OPTION__PHOTOELECTRONSPECTRUM__PES_FLUX)
      call pes_flux_map_from_states(pflux, restart, st, llg, pesP, krng, Lp)      
    end select

    call output_spin_pes()
    
  end if



  write(message(1), '(a)') 'Done'
  call messages_info(1)

  call messages_print_stress(stdout)

  call restart_end(restart)    

  call states_end(st)

  call geometry_end(geo)
  call simul_box_end(sb)
  call space_end(space)

  call io_end()
  call messages_end()

  call parser_end()
  call global_end()
  
  SAFE_DEALLOCATE_A(pesP)    
  SAFE_DEALLOCATE_A(pmesh)
  SAFE_DEALLOCATE_P(Lp)
  SAFE_DEALLOCATE_A(Ekin)
  if (.not. need_pmesh .or. pes_method == OPTION__PHOTOELECTRONSPECTRUM__PES_MASK) then
    SAFE_DEALLOCATE_P(Lg)
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
        SAFE_ALLOCATE(pesP_out(1:llp(1),1:llp(2),1:llp(3)))
        pesP_out(:,:,:) = pesP(:,:,:,1)+pesP(:,:,:,2)
    
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
      
      if(bitand(pesout%what, OPTION__PHOTOELECTRONSPECTRUMOUTPUT__ENERGY_TOT) /= 0) then
        call messages_print_stress(stdout, "Energy-resolved PES")

        select case (pes_method)
        case (OPTION__PHOTOELECTRONSPECTRUM__PES_MASK)
          call pes_mask_output_power_totalM(pesP_out,outfile('./PES_power',ist, ispin, 'sum'), &
                                            Lg, llp, dim, Emax, Estep, interpolate = .true.)
        case (OPTION__PHOTOELECTRONSPECTRUM__PES_FLUX)                                     
          call messages_not_implemented("Energy-resolved PES for the flux method") 
        end select 
        
      end if
      
      if(bitand(pesout%what, OPTION__PHOTOELECTRONSPECTRUMOUTPUT__ENERGY_ANGLE) /= 0) then
        call messages_print_stress(stdout, "Angle- and energy-resolved PES")
        
        select case (pes_method)
        case (OPTION__PHOTOELECTRONSPECTRUM__PES_MASK)
          call pes_mask_output_ar_polar_M(pesP_out,outfile('./PES_angle_energy',ist, ispin, 'map'), &
                                          Lg, llp, dim, pol, Emax, Estep)
        case (OPTION__PHOTOELECTRONSPECTRUM__PES_FLUX)                                     
          call messages_not_implemented("Angle- and energy-resolved PES for the flux method") 
        end select                       
                                          
                                          
      end if

      if(bitand(pesout%what, OPTION__PHOTOELECTRONSPECTRUMOUTPUT__VELOCITY_MAP_CUT) /= 0) then
        call messages_print_stress(stdout, "Velocity map on a plane")
        dir = -1
        if(sum((pvec-(/1 ,0 ,0/))**2)  <= M_EPSILON  )  dir = 1
        if(sum((pvec-(/0 ,1 ,0/))**2)  <= M_EPSILON  )  dir = 2
        if(sum((pvec-(/0 ,0 ,1/))**2)  <= M_EPSILON  )  dir = 3

        if (have_zweight_path) then
          filename = outfile('PES_velocity_map',ist,ispin,'path')
        else
          filename = outfile('PES_velocity_map',ist,ispin,'p'//index2axis(dir)//'=0')
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
          call pes_out_velocity_map_cut(pesP_out, filename, llp, dim, pol, dir, integrate, &
                                             pos = idxZero, pmesh = pmesh)
!           call pes_mask_output_full_mapM_cut(pesP_out, filename, llp, dim, pol, dir, integrate, &
!                                              pos = idxZero, pmesh = pmesh)
        else
          call pes_out_velocity_map_cut(pesP_out, filename, llp, dim, pol, dir, integrate, &
                                             pos = idxZero, Lk = Lg)
!           call pes_mask_output_full_mapM_cut(pesP_out, filename, llp, dim, pol, dir, integrate, &
!                                              pos = idxZero, Lk = Lg)
        end if
      end if

      if(bitand(pesout%what, OPTION__PHOTOELECTRONSPECTRUMOUTPUT__ENERGY_XY) /= 0) then
        call messages_print_stress(stdout, "Angle and energy-resolved on a plane")
        if(uEstep >  0 .and. uEstep > Estep) then
          Estep = uEstep
        else
          Estep = Emax/size(Lg,1)
        end if

        select case (pes_method)
        case (OPTION__PHOTOELECTRONSPECTRUM__PES_MASK)
          call pes_mask_output_ar_plane_M(pesP_out,outfile('./PES_energy',ist,ispin,'map'), &
                                          Lg, llp, dim, pol, Emax, Estep)
        case (OPTION__PHOTOELECTRONSPECTRUM__PES_FLUX)                                     
          call messages_not_implemented("Angle and energy-resolved on a plane for the flux method") 
        end select   
                                        
      end if

      if(bitand(pesout%what, OPTION__PHOTOELECTRONSPECTRUMOUTPUT__ENERGY_TH_PH) /= 0) then
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
          call pes_mask_output_ar_spherical_cut_M(pesP_out,outfile('./PES_sphere',ist,ispin,'map'), & 
                                                  Lg, llp, dim, pol, Emin, Emax, Estep)

        case (OPTION__PHOTOELECTRONSPECTRUM__PES_FLUX)                                     
          call messages_not_implemented("PES on spherical cuts for the flux method") 
        end select                                          

      end if

      if(bitand(pesout%what, OPTION__PHOTOELECTRONSPECTRUMOUTPUT__VELOCITY_MAP) /= 0) then
        
        call io_function_read_how(sb, namespace, how, ignore_error = .true.)
        call messages_print_stress(stdout, "Full velocity map")
        
        filename = outfile('./PES_velocity_map', ist, ispin)
        if (need_pmesh) then
          !force vtk output
          how = io_function_fill_how("VTK")
           
!           call pes_mask_output_full_mapM(pesP_out, filename, Lg, llp, how, sb, pmesh)
          call pes_out_velocity_map(pesP_out, filename, Lg, llp, how, sb, pmesh)
        else
!           call pes_mask_output_full_mapM(pesP_out, filename, Lg, llp, how, sb)
          call pes_out_velocity_map(pesP_out, filename, Lg, llp, how, sb)
        end if
        
      end if

      if(bitand(pesout%what, OPTION__PHOTOELECTRONSPECTRUMOUTPUT__ARPES) /= 0) then
        call messages_print_stress(stdout, "ARPES")

        forall (i1=1:llp(1), i2=1:llp(2), i3=1:llp(3))
          pmesh(i1,i2,i3,dim) = units_from_atomic(units_out%energy, &
            sign(M_ONE,pmesh(i1,i2,i3,dim)) * sum( pmesh(i1,i2,i3,1:dim)**2 )/M_TWO)
        end forall

        how = io_function_fill_how("VTK")

!         call pes_mask_output_full_mapM(pesP_out, outfile('./PES_ARPES', ist, ispin), &
!                                        Lg, llp, how, sb, pmesh)
        call pes_out_velocity_map(pesP_out, outfile('./PES_ARPES', ist, ispin), &
                                       Lg, llp, how, sb, pmesh)
      end if
      
      
      if(bitand(pesout%what, OPTION__PHOTOELECTRONSPECTRUMOUTPUT__ARPES_CUT) /= 0) then
        call messages_print_stress(stdout, "ARPES cut on reciprocal space path")
        
        filename = outfile('./PES_ARPES', ist, ispin, "path")
        call pes_out_arpes_cut(pesP_out, filename, llp, pmesh, Ekin)
        
      end if
      
      
    end subroutine output_pes
    
    subroutine get_kpath_direction(kpoints, krng, kpth_dir, pvec)
      type(kpoints_t),   intent(in)  :: kpoints 
      integer,           intent(in)  :: krng(2)
      integer,           intent(out) :: kpth_dir
      FLOAT,             intent(out) :: pvec(3)

         
      FLOAT                :: kpt(3)
         
      PUSH_SUB(get_kpath_direction)
      
      kpth_dir = -1
         
      kpt = M_ZERO
      kpt(1:dim) = kpoints_get_point(kpoints, krng(1)+1)-kpoints_get_point(kpoints, krng(1))
      kpt(1:dim) = kpt(1:dim)/sqrt(sum(kpt(1:dim)**2))  
           
      
      if (sum((kpt(:) - (/1,0,0/))**2) < M_EPSILON) then
        kpth_dir = 1
        write(message(1), '(a)') 'along kx'
        pvec = (/0,1,0/)
      end if        
              
      if (sum((kpt(:) - (/0,1,0/))**2) < M_EPSILON) then 
        kpth_dir = 2
        write(message(1), '(a)') 'along ky'
        pvec = (/1,0,0/)        
      end if
      
      call messages_info(1)
      
    
      if (kpth_dir == -1) then
        message(1) = "K-points with zero weight path works only with paths along kx or ky."
        call messages_fatal(1)
      end if
      
      POP_SUB(get_kpath_direction)      
    end subroutine get_kpath_direction

    
    subroutine write_kpoints_info(kpoints, ikstart, ikend)
      type(kpoints_t),   intent(in) :: kpoints 
      integer,           intent(in) :: ikstart   
      integer,           intent(out):: ikend  

      integer :: ik, idir
      character(len=100) :: str_tmp
            
      PUSH_SUB(write_kpoints_info)

      do ik = ikstart, ikend
        write(message(1),'(i8,1x)') ik
        write(str_tmp,'(f12.6)') kpoints_get_weight(kpoints, ik)
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
        
        no_l = 0
        if(parse_block(namespace, 'TDExternalFields', blk) == 0) then
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
