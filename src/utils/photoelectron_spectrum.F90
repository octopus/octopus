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
!! $Id$

#include "global.h"

program photoelectron_spectrum
  use command_line_m
  use geometry_m
  use global_m
  use grid_m
  use kpoints_m
  use io_binary_m
  use io_function_m
  use io_m
  use loct_m
  use messages_m
  use parser_m
  use pes_m  
  use pes_mask_m  
  use profiling_m
  use restart_m
  use simul_box_m
  use sort_om
  use space_m
  use string_m
  use states_m
  use unit_m
  use unit_system_m
  use utils_m
  use varinfo_m
  
  implicit none

  integer              :: ierr, mode, interp, integrate
  integer              :: dim, ll(3), lll(3), dir, how, idim
  FLOAT                :: Emax, Emin, Estep, uEstep,uEspan(2), pol(3)
  FLOAT                :: uThstep, uThspan(2), uPhstep, uPhspan(2), pvec(3)
  FLOAT                :: center(3)
  FLOAT, pointer       :: Lk(:,:), RR(:)
  FLOAT, allocatable   :: pesk(:,:,:), pmesh(:,:,:,:)
  integer, allocatable :: Lp(:,:,:,:,:)
  logical              :: interpol, need_pmesh
  integer              :: ii, i1,i2,i3, idxZero(1:3)
  
  type(space_t)     :: space
  type(geometry_t)  :: geo
  type(simul_box_t) :: sb
  type(states_t)    :: st
  type(grid_t)      :: gr
  type(restart_t)   :: restart
  
  character(len=512) :: filename


  call getopt_init(ierr)
  if(ierr /= 0) then
    message(1) = "Your Fortran compiler doesn't support command-line arguments;"
    message(2) = "the oct-photoelectron-spectrum command is not available."
    call messages_fatal(2)
  end if
!   ! This first time checks only  if the --help option is present
!   call getopt_photoelectron_spectrum(mode,interp,uEstep, uEspan,&
!                                      uThstep, uThspan, uPhstep, &
!                                      uPhspan, pol, center, pvec, integrate)
!


  call global_init(is_serial = .true.)
  
  call messages_init()  
  call io_init()

  !* In order to initialize k-points
  call unit_system_init()
  
  call space_init(space)
  call geometry_init(geo, space)
  call simul_box_init(sb, geo, space)
  gr%sb = sb
  call states_init(st, gr, geo)
  !*


  !Initial values
  ll(:) = 1 
  interpol = .true. 
  need_pmesh = (kpoints_number(sb%kpoints) > 1)

  !set default values
  mode = 1
  interp = 1
  integrate = -1
  uEstep = -1
  uEspan = (/-1,-1/)
  uThstep = -1
  uThspan = (/-1,-1/)
  uPhstep = -1
  uPhspan = (/-1,-1/)
  center = (/0,0,0/)
  pvec = (/1,0,0/)
  Emin = M_ZERO
  Emax = M_ZERO
  
  if(simul_box_is_periodic(sb)) then
    mode = 7
!     if(sb%dim == 2 ) mode = 3
!     if(sb%dim == 3 ) mode = 7
  end if

  
  call get_laser_polarization(pol)
  if (sum(pol(1:3)**2) <= M_EPSILON) pol = (/0,0,1/)
  
  
  call getopt_photoelectron_spectrum(mode,interp,uEstep, uEspan,&
                                     uThstep, uThspan, uPhstep, &
                                     uPhspan, pol, center, pvec, integrate)
                                     
  call getopt_end()
                                       
  if(interp  ==  0) interpol = .false.

  call messages_print_stress(stdout)
  call pes_mask_read_info("td.general/", dim, Emax, Estep, ll(:), Lk, RR)


  write(message(1), '(a)') 'Read PES info file.'
  call messages_info(1)
  
  
  lll(:) = ll(:)
  call restart_module_init()
  call restart_init(restart, RESTART_TD, RESTART_TYPE_LOAD, st%dom_st_kpt_mpi_grp, ierr)
  if(ierr /= 0) then
    message(1) = "Unable to read time-dependent restart information."
    call messages_fatal(1)
  end if
  
  SAFE_ALLOCATE(Lp(1:ll(1),1:ll(2),1:ll(3),kpoints_number(sb%kpoints),1:3))

  ll(1:sb%dim) = ll(1:sb%dim) * sb%kpoints%nik_axis(1:sb%dim)

  SAFE_ALLOCATE(pmesh(1:ll(1),1:ll(2),1:ll(3),1:3))
  call pes_mask_pmesh(sb%kpoints, lll, Lk, pmesh, Lp, idxZero)  
 
  SAFE_ALLOCATE(pesk(1:ll(1),1:ll(2),1:ll(3)))
  
!   do i1 = 1, ll(1)
!     do i2 = 1, ll(2)
!       do i3 = 1, ll(3)
!         if ( (pmesh(i1, i2, i3, 1) - M_HUGE)**2<=M_EPSILON ) then
!           print *,"unfilled indices (ip1,ip2,ip3)=     ", i1,i2,i3
!         end if
!
!       end do
!     end do
!   end do


  call pes_mask_map_from_states(restart, st, lll, pesk, Lp)

  call restart_end(restart)    
  
  if (.not. simul_box_is_periodic(sb) .or. kpoints_number(sb%kpoints) == 1) then
    ! There is no need to use pmesh we just need to sort Lk in place
    ! in order to have a coordinate ordering coherent with pesk
    do idim = 1, sb%dim
      call sort(Lk(1:ll(idim), idim)) 
    end do  
  end if  
    
!   if(.false.)
!     !Read directly from the obf file
!     filename=io_workpath('td.general/PESM_map.obf')
!     call io_binary_read(trim(filename),ll(1)*ll(2)*ll(3),pesk, ierr)
!     if(ierr > 0) then
!       message(1) = "Failed to read file "//trim(filename)
!       call messages_fatal(1)
!     end if
!   end if


  write(message(1), '(a)') 'Read PES restart files.'
  call messages_info(1)

  !! set user values
  if(uEstep >  0 .and. uEstep > Estep)    Estep = uEstep
  if(uEspan(1) > 0 ) Emin = uEspan(1)
  if(uEspan(2) > 0 ) Emax = uEspan(2)


  call unit_system_init()
  call messages_print_stress(stdout)
 
  write(message(1),'(a,f10.2,a2,f10.2,a2,f10.2,a1)') &
                   "Zenith axis: (",pol(1),", ",pol(2),", ",pol(3),")"
  call messages_info(1)


  ! Convert the grid units
  if (need_pmesh) then
    
    forall (i1=1:ll(1), i2=1:ll(2), i3=1:ll(3), ii = 1:3)
      pmesh(i1,i2,i3,ii) = units_from_atomic(sqrt(units_out%energy), pmesh(i1,i2,i3,ii))
    end forall
  end if

  ! choose what to calculate
  ! these functions are defined in pes_mask_out_inc.F90
  select case(mode)
  case(1) ! Energy-resolved
    write(message(1), '(a)') 'Save energy-resolved PES'
    call messages_info(1)
    call pes_mask_output_power_totalM(pesk,'./PES_power.sum', Lk, ll, dim, Emax, Estep, interpol)
 
 
  case(2) ! Angle and energy resolved
    write(message(1), '(a)') 'Save angle- and energy-resolved PES'
    call messages_info(1)
    call pes_mask_output_ar_polar_M(pesk,'./PES_angle_energy.map', Lk, ll, dim, pol, Emax, Estep)


  case(3) ! On a plane
    
    dir = -1
    if(sum((pvec-(/1 ,0 ,0/))**2)  <= M_EPSILON  )  dir = 1
    if(sum((pvec-(/0 ,1 ,0/))**2)  <= M_EPSILON  )  dir = 2
    if(sum((pvec-(/0 ,0 ,1/))**2)  <= M_EPSILON  )  dir = 3

    filename = "PES_velocity_map.p"//index2axis(dir)//"=0"


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
      filename = "PES_velocity.map.i_"//trim(index2var(integrate))//".p"//index2axis(dir)//"=0"
    end if
    
    if (need_pmesh) then
      call pes_mask_output_full_mapM_cut(pesk, filename, ll, dim, pol, dir, integrate, & 
                                         pos = idxZero, pmesh = pmesh)  
    else    
      call pes_mask_output_full_mapM_cut(pesk, filename, ll, dim, pol, dir, integrate, &
                                         pos = idxZero, Lk = Lk)    
    end if

  case(4) ! Angle energy resolved on plane 
    write(message(1), '(a)') 'Save angle and energy-resolved PES'
    call messages_info(1)
    if(uEstep >  0 .and. uEstep > Estep) then
      Estep = uEstep
    else
      Estep = Emax/size(Lk,1)
    end if

    call pes_mask_output_ar_plane_M(pesk,'./PES_energy.map', Lk, ll, dim, pol, Emax, Estep)

  case(5) ! Angular-resolved  


    write(message(1), '(a,es19.12,a2,es19.12,2x,a19)') &
          'Save PES on a spherical cut at E= ',Emin,", ",Emax, & 
           str_center('['//trim(units_abbrev(units_out%energy)) // ']', 19) 
    call messages_info(1)

    if(uEstep >  0 .and. uEstep > Estep) then
     Estep = uEstep
    else
     Estep = Emax/size(Lk,1)
    end if
 
    call pes_mask_output_ar_spherical_cut_M(pesk,'./PES_sphere.map', Lk, ll, dim, pol, Emin, Emax, Estep)

  case(6) ! Full momentum resolved matrix 
  
    call io_function_read_how(sb, how, ignore_error = .true.)
 
    write(message(1), '(a)') 'Save full momentum-resolved PES'
    call messages_info(1)

    if (need_pmesh) then
      how = io_function_fill_how("VTK")
      call pes_mask_output_full_mapM(pesk, './PES_velocity_map', Lk, ll, how, sb, pmesh)
    else
      call pes_mask_output_full_mapM(pesk, './PES_velocity_map', Lk, ll, how, sb) 
    end if

    case(7) ! ARPES 
 
      write(message(1), '(a)') 'Save ARPES'
      call messages_info(1)

      forall (i1=1:ll(1), i2=1:ll(2), i3=1:ll(3))
        pmesh(i1,i2,i3,sb%dim) = units_from_atomic(units_out%energy, &
          sign(M_ONE,pmesh(i1,i2,i3,sb%dim)) * sum( pmesh(i1,i2,i3,1:sb%dim)**2 )/M_TWO)
      end forall

      how = io_function_fill_how("VTK")
      
      call pes_mask_output_full_mapM(pesk, './PES_ARPES', Lk, ll, how, sb, pmesh)   

  end select


  write(message(1), '(a)') 'Done'
  call messages_info(1)

  call messages_print_stress(stdout)

  call states_end(st)

  call geometry_end(geo)
  call simul_box_end(sb)
  call space_end(space)

  call io_end()
  call messages_end()
  call global_end()
  
  SAFE_DEALLOCATE_A(pesk)    
  SAFE_DEALLOCATE_A(pmesh)
  SAFE_DEALLOCATE_A(Lp)
  SAFE_DEALLOCATE_P(Lk)
  
  contains

    subroutine get_laser_polarization(lPol)
       FLOAT,   intent(out) :: lPol(:) 
       
        type(block_t)       :: blk
        integer             :: no_l
        
        PUSH_SUB(get_laser_polarization)
        
        no_l = 0
        if(parse_block('TDExternalFields', blk) == 0) then
          no_l = parse_block_n(blk)

          call parse_block_float(blk, 0, 1, lPol(1))
          call parse_block_float(blk, 0, 2, lPol(2))
          call parse_block_float(blk, 0, 3, lPol(3))


          call parse_block_end(blk)
        end if
        
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

    !< Generate the momentum-space mesh (p) and the arrays mapping the 
    !< the mask and the kpoint meshes in p.
    subroutine pes_mask_pmesh(kpoints, ll, Lk, pmesh, Lp, idxZero)
      type(kpoints_t),   intent(inout) :: kpoints 
      integer,           intent(in)    :: ll(:)             !< ll(1:dim): the dimensions of the mask-mesh
      FLOAT,             intent(in)    :: Lk(:,:)           !< Lk(1:maxval(ll),1:dim): the mask-mesh points  
      FLOAT,             intent(out)   :: pmesh(:,:,:,:)    !< pmesh(i1,i2,i3,1:dim): contains the positions of point 
                                                          !< in the final mesh in momentum space "p" combining the 
                                                          !< mask-mesh with kpoints. 
      integer,           intent(out) :: Lp(:,:,:,:,:)     !< Lp(1:ll(1),1:ll(2),1:ll(3),1:nkpt,1:dim): maps a 
                                                          !< mask-mesh triplet of indices together with a kpoint 
                                                          !< index into a triplet on the combined momentum space mesh.

      integer,           intent(out) :: idxZero(:)       !< The triplet identifying the zero of the coordinates                  
  
      integer :: ik, i1, i2, i3, nk(1:3), ip1, ip2, ip3, idir, dim
      FLOAT :: kpt(1:3),kval(1:3), dx(1:sb%dim)
      integer, allocatable :: Lkpt(:,:), err
      
      integer, allocatable :: idx(:,:)
      FLOAT, allocatable   :: Lk_(:,:)
      
  
      PUSH_SUB(pes_mask_pmesh)

      SAFE_ALLOCATE(Lkpt(1:kpoints_number(kpoints),1:3))
      
      dim = sb%dim
      
      nk(:) = 1  
      nk(1:dim) = kpoints%nik_axis(1:dim)

      Lkpt(:,:) = 1
      kpt(:) = M_ZERO
            
      call kpoints_grid_generate(dim, kpoints%nik_axis(1:dim), kpoints%shifts(1:dim), &
                                 kpoints%full%red_point,  Lkpt(:,1:dim))

!       do ik = 1, kpoints_number(kpoints)
!         kpt(1:sb%dim) = kpoints_get_point(kpoints, ik, absolute_coordinates = .false.)
!         print *, ik, "Lkpt(ik)= ", Lkpt(ik,:), "-- kpt= ",kpt
!       end do
      
      
!       print *,"----"
!       print *,"ll(:)", ll(:)
!       print *,"----"
      
      
      ! We want the results to be sorted on a cube i,j,k
      ! with the first triplet associated with the smallest positions 
      SAFE_ALLOCATE(idx(1:maxval(ll(:)), 1:3))
      SAFE_ALLOCATE(Lk_(1:maxval(ll(:)), 1:3))
      idx(:,:)=1
      do idir = 1, dim
        Lk_(:,idir) = Lk(:,idir)
        call sort(Lk_(1:ll(idir), idir), idx(1:ll(idir), idir)) 
      end do  
      
      pmesh(:, :, :, :) = M_HUGE      
      err = -1
      
      ! Generate the p-space mesh and populate Lp
      do ik = 1, kpoints_number(kpoints)
        kpt(1:dim) = kpoints_get_point(kpoints, ik) 
        do i1 = 1, ll(1) 
          do i2 = 1, ll(2) 
            do i3 = 1, ll(3) 
              
              kval (1:3)= (/Lk_(i1,1),Lk_(i2,2),Lk_(i3,3)/)               

              ip1 = (i1 - 1) * nk(1) + Lkpt(ik,1) 
              ip2 = (i2 - 1) * nk(2) + Lkpt(ik,2) 
              ip3 = (i3 - 1) * nk(3) + Lkpt(ik,3)
              
              Lp(idx(i1,1),idx(i2,2),idx(i3,3),ik,1:3) =  (/ip1,ip2,ip3/)
              
              
              pmesh(ip1, ip2, ip3, 1:dim) = kval(1:dim) + kpt(1:dim)

!               print *,ik,i1,i2,i3,"  Lp(i1,i2,i3,ik,1:dim) = ",  (/ip1,ip2,ip3/)
!               print *, "pmesh = ",pmesh(ip1, ip2, ip3, 1:dim) ,"  kval = ",  kval (1:3), "  kpt = ", kpt(1:dim)

              if (sum(pmesh(ip1, ip2, ip3, 1:3)**2)<=M_EPSILON) then
                err = err + 1 
                idxZero(1:3) = (/ip1,ip2,ip3/)
              end if


            end do 
          end do 
        end do 
      end do
  
      if (err == -1) then
        write(message(1), '(a)') 'Illformed momentum-space mesh: could not find zero coordinate.'
        call messages_fatal(1)
      end if 

      if (err > 0) then
        write(message(1), '(a)') 'Illformed momentum-space mesh: too many points with zero coordinate .'
        call messages_fatal(1)
      end if 

  
  
      SAFE_DEALLOCATE_A(Lkpt)
      
      POP_SUB(pes_mask_pmesh)
    end subroutine pes_mask_pmesh


    !< Build the photoemission map form the restart files
    subroutine pes_mask_map_from_states(restart, st, ll, pesK, Lp)
      type(restart_t),    intent(in) :: restart
      type(states_t),     intent(in) :: st
      integer,            intent(in) :: ll(:)
      FLOAT, target,     intent(out) :: pesK(:,:,:)
      integer, optional,  intent(in) :: Lp(:,:,:,:,:)
      
      integer :: ik, ist, idim, itot
      integer :: i1, i2, i3, ip(1:3)
      integer :: idone, ntodo
      CMPLX   :: psik(1:ll(1),1:ll(2),1:ll(3))

      PUSH_SUB(pes_mask_map_from_states)
  
      ntodo = st%d%kpt%nglobal * st%nst * st%d%dim
      idone = 0 
      call loct_progress_bar(-1, ntodo)
      
      pesK = M_ZERO
      do ik = 1, st%d%kpt%nglobal
        do ist = 1, st%nst
          do idim = 1, st%d%dim
            itot = ist + (ik-1) * st%nst +  (idim-1) * st%nst * st%d%kpt%nglobal
            call pes_mask_map_from_state(restart, itot, ll, psik)
            
            do i1=1, ll(1)
              do i2=1, ll(2)
                do i3=1, ll(3)
                  ip(1:3) = Lp(i1 , i2, i3, ik, 1:3) 
                  
                  pesK(ip(1),ip(2),ip(3)) = pesK(ip(1),ip(2),ip(3)) &
                    + abs(psik(i1,i2,i3))**2 * st%occ(ist, ik) * st%d%kweights(ik)
                                          
                end do
              end do
            end do
            
            idone = idone +1 
            call loct_progress_bar(idone, ntodo)
            
          end do
        end do
      end do

      write(stdout, '(1x)')

      POP_SUB(pes_mask_map_from_states)
    end subroutine pes_mask_map_from_states


    subroutine pes_mask_map_from_state(restart, idx, ll, psik)
      type(restart_t),  intent(in)  :: restart
      integer,          intent(in)  :: idx
      integer,          intent(in)  :: ll(:)
      CMPLX, target,    intent(out) :: psik(:,:,:)

      character(len=80) :: filename, path
      integer ::  np, err, iunit 
      character(len=128) :: lines(2)
  
      PUSH_SUB(pes_mask_map_from_state)

      psik = M_Z0
      np = product(ll(:))
      
      write(filename,'(i10.10)') idx
 
      path = trim(restart_dir(restart))//'/pes_'//trim(filename)//'.obf'
      call io_binary_read(path, np, psik(:,:,:), err)
      if (err /= 0) then
        message(1) = "Unable to read PES mask restart data from '"//trim(path)//"'."
        call messages_warning(1)
      end if

      POP_SUB(pes_mask_map_from_state)  
    end subroutine pes_mask_map_from_state
    



end program photoelectron_spectrum

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
