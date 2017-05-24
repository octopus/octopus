!! Copyright (C) 2017 U. De Giovannini
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

program floquet_observables
  use calc_mode_par_oct_m
  use comm_oct_m
  use command_line_oct_m
  use geometry_oct_m
  use fft_oct_m
  use floquet_oct_m
  use global_oct_m
  use grid_oct_m
  use hamiltonian_oct_m
  use kpoints_oct_m
  use io_binary_oct_m
  use io_function_oct_m
  use io_oct_m
  use loct_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use messages_oct_m
  use mpi_oct_m
  use multicomm_oct_m
  use parser_oct_m
  use profiling_oct_m
  use restart_oct_m
  use simul_box_oct_m
  use system_oct_m
  use sort_oct_m
  use space_oct_m
  use string_oct_m
  use states_oct_m
  use states_dim_oct_m
  use states_io_oct_m
  use types_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use utils_oct_m
  use varinfo_oct_m
  
  implicit none

  integer              :: out_what  

  integer              :: ierr
  integer              :: dim, dir, how, idim, pdim
  integer              :: ii, i1,i2,i3
  type(block_t)        :: blk  
  
  type(restart_t)      :: restart
  type(system_t)      :: sys
  type(hamiltonian_t) :: hm
  character(len=512)   :: filename, str, str2
  integer              :: ist, ispin  
  type(states_t)          :: dressed_st 
  type(states_t), pointer ::    gs_st

  call getopt_init(ierr)
  if(ierr /= 0) then
    message(1) = "Your Fortran compiler doesn't support command-line arguments;"
    message(2) = "the oct-photoelectron-spectrum command is not available."
    call messages_fatal(2)
  end if


  call global_init()
  
  call messages_init()  
  call io_init()
  call calc_mode_par_init()

  call fft_all_init()
  call unit_system_init()
  
  call system_init(sys)
  

  call floquet_init(sys,hm%F,sys%st%d%dim)
  gs_st => sys%st

  call states_init(dressed_st, sys%gr, sys%geo,floquet_dim=hm%F%floquet_dim)
  call kpoints_distribute(dressed_st%d,sys%mc)
  call states_distribute_nodes(dressed_st,sys%mc)
  call states_exec_init(dressed_st, sys%mc)
  call restart_module_init()

  !%Variable FloquetObservableCalc
  !%Type flag
  !%Default none
  !%Section Utilities::oct-floquet_observables
  !%Description
  !% Specifies what observables to calculate form the Floquet solution.
  !% Example: <tt>f_norms + f_spin</tt>
  !%Option f_norms bit(1)
  !% The norms of the Floquet states in each tensorial subspace.
  !%Option f_arpes bit(2)
  !% Calculate ARPES matrix elements for Floquet states.
  !%Option f_spin bit(3)
  !% Calculate the spin polarization of each state. (Not implemented)
  !%Option f_td_spin bit(4)
  !% Calculate the time-dependent spin projections of the Floquet eigenstates
  !% given by FloquetObservableTDspinKpoints and FloquetObservableTDspinState
  !%Option f_hhg bit(5)
  !% Calculate the HHG spectrum.
  !%End
  call parse_variable('FloquetObservableCalc', out_what, out_what)
  
!   call getopt_floquet_observables(uEstep, uEspan,&
!                                      uThstep, uThspan, uPhstep, &
!                                      uPhspan, pol, center, pvec, integrate)
!
!   call getopt_end()
                                       

  ! load floquet states only for certain tasks
  if(.not. (iand(out_what, OPTION__FLOQUETOBSERVABLECALC__F_TD_SPIN) /= 0) ) then
       call states_allocate_wfns(dressed_st,sys%gr%der%mesh, wfs_type = TYPE_CMPLX)
       call floquet_restart_dressed_st(hm, sys, dressed_st, ierr)
       call messages_write('Read Floquet restart files.')
       call messages_info()
  end if

  if(iand(out_what, OPTION__FLOQUETOBSERVABLECALC__F_NORMS) /= 0) then
    call messages_write('Calculate norms of Floquet subspaces.')
    call messages_info()

    call calc_floquet_norms(sys%gr%der%mesh,sys%gr%sb%kpoints,gs_st,dressed_st, hm%F%iter,hm%F%floquet_dim)
  end if

  if(iand(out_what, OPTION__FLOQUETOBSERVABLECALC__F_ARPES) /= 0) then
    call messages_write('Calculate Floqeut ARPES.')
    call messages_info()

    call calc_floquet_arpes()
  end if

  if(iand(out_what, OPTION__FLOQUETOBSERVABLECALC__F_TD_SPIN) /= 0) then
    call messages_write('Calculate td-spin of Floquet eigenstates.')
    call messages_info()

    if(gs_st%d%dim /=2 ) then
      call messages_write('Need spin resolved calculation')
      call messages_fatal()
    end if

    call calc_floquet_td_spin()
  end if
  
  

  call system_end(sys)
  call fft_all_end()
  call io_end()
  call messages_end()
  call global_end()

contains 
    
  subroutine calc_floquet_arpes()
  
  FLOAT, allocatable :: spect(:,:), me(:,:) 
  FLOAT :: pomega, pol(1:3)
  type(block_t)        :: blk  
  
  
  !%Variable FloquetObservablePesOmega
  !%Type float
  !%Default 50 eV
  !%Section Utilities::oct-floquet_observables
  !%Description
  !% The probe energy needed to calculate  photoemission matrix elements.
  !%End
  call parse_variable('FloquetObservablePesOmega', CNST(1.83749219065), pomega, units_inp%energy)
  call messages_print_var_value(stdout,'Frequency of PES probe field', pomega)
  
  
  !%Variable FloquetObservablePesPol
  !%Type block
  !%Default 
  !%Section Utilities::oct-floquet_observables
  !%Description
  !% Probe field polarization direction.
  !%
  !% <tt>%FloquetObservablePesPol
  !% <br> px | py | pz
  !% <br>%</tt>
  !%End
  pol(:)=M_ZERO
  pol(1)=M_ONE

  if(parse_block('FloquetObservablePesPol', blk) == 0) then
    if(parse_block_cols(blk,0) /= sys%gr%sb%dim) call messages_input_error('FloquetObservablePesPol')
    do idim = 1, sys%gr%sb%dim
      call parse_block_float(blk, 0, idim - 1, pol(idim))
    end do
    call parse_block_end(blk)
  end if
  
  write(message(1),'(a,f4.2,a,f4.2,a,f4.2,a)') 'Info: ARPES probe polarization: (', pol(1),',', pol(2),',', pol(3),')'
  call messages_info(1)
    
    
  SAFE_ALLOCATE(spect(dressed_st%nst,dressed_st%d%nik))
  SAFE_ALLOCATE(me(dressed_st%nst,dressed_st%d%nik))
  call floquet_photoelectron_spectrum(hm, sys, dressed_st, pomega, pol, spect, me)

  write(str,'(I5)') hm%F%iter
  if(simul_box_is_periodic(sys%gr%sb)) then
    filename = 'floquet_arpes_me_'//trim(adjustl(str))
    call states_write_bandstructure(FLOQUET_DIR, dressed_st%nst, dressed_st, sys%gr%sb, filename, vec = me)

    filename = 'floquet_arpes_'//trim(adjustl(str))
    call states_write_bandstructure(FLOQUET_DIR, dressed_st%nst, dressed_st, sys%gr%sb, filename, vec = spect)
  end if

  SAFE_DEALLOCATE_A(spect)
  SAFE_DEALLOCATE_A(me)
    
  end subroutine calc_floquet_arpes  

  !-------------------------------------------------
  subroutine calc_floquet_td_spin()

  FLOAT :: spin(1:3)
  type(block_t)        :: blk
  integer :: nk_input, ik,iik, ist,iist, nst_input, it, iunit
  integer :: nik, dim, nst, itot
  integer, allocatable :: kpoints_input(:), states_input(:)
  CMPLX, allocatable   :: psi_t(:,:), zpsi(:), F_psi(:,:)
  FLOAT :: time

  !%Variable FloquetObservableTDspinKpoints
  !%Type block
  !%Default none
  !%Section Utilities::oct-floquet_observables
  !%Description
  !% Indices of k-points for td-spin calculation  
  !%End

  if(parse_block('FloquetObservableTDspinKpoints', blk) == 0) then
    nk_input = parse_block_cols(blk, 0)
    if(nk_input < 1 ) call messages_input_error('FloquetObservableTDspinKpoints')
    SAFE_ALLOCATE(kpoints_input(nk_input))
    do ik = 1, nk_input
      call parse_block_integer(blk, 0, ik-1, kpoints_input(ik))
    end do
    call parse_block_end(blk)
  end if
  write(message(1),'(a,i3,a)') 'Info: Read ', nk_input ,' kpoints for td-spin calculation'
  call messages_info(1)
  
  !%Variable FloquetObservableTDspinStates
  !%Type block
  !%Default none
  !%Section Utilities::oct-floquet_observables
  !%Description
  !% Indices of Floquet eigenstates for td-spin calculation
  !%End

  if(parse_block('FloquetObservableTDspinStates', blk) == 0) then
    nst_input = parse_block_cols(blk, 0)
    if(nst_input < 1 ) call messages_input_error('FloquetObservableTDspinStates')
    SAFE_ALLOCATE(states_input(nst_input))
    do ist = 1, nst_input
      call parse_block_integer(blk, 0, ist-1, states_input(ist))
    end do
    call parse_block_end(blk)
  end if
  write(message(1),'(a,i3,a)') 'Info: Read ', nst_input ,' states for td-spin calculation'
  call messages_info(1)

  ! prepare restart structure
  call restart_init(restart, RESTART_FLOQUET, RESTART_TYPE_LOAD, &
                       sys%mc, ierr,sys%gr%der%mesh)
  call states_look(restart, nik, dim, nst, ierr)
  if(dim/=dressed_st%d%dim .or. nik/=sys%gr%sb%kpoints%reduced%npoints.or.nst/=dressed_st%nst) then
     write(message(1),'(a)') 'Did not find commensurate Floquet restart structure'
     call messages_fatal(1)
  end if

  SAFE_ALLOCATE(psi_t(1:sys%gr%der%mesh%np,sys%st%d%dim))
  SAFE_ALLOCATE(zpsi(1:sys%gr%mesh%np))
  SAFE_ALLOCATE(F_psi(1:sys%gr%mesh%np,dressed_st%d%dim))

  do ik=1,nk_input
    iik=kpoints_input(ik)
    write(str,'(I5)') iik
    
    do ist=1,nst_input
      iist=states_input(ist)
      write(str2,'(I5)') iist

      ! read the floquet wavefunction for states iik,iist
      do idim = 1, dressed_st%d%dim
         itot = idim + (iist-1)*dressed_st%d%dim +  (iik-1)*dressed_st%nst*dressed_st%d%dim 
         write(filename,'(i10.10)') itot
         call zrestart_read_mesh_function(restart, trim(adjustl(filename)), sys%gr%mesh, zpsi, ierr)
         print *, ierr, trim(adjustl(filename)), iik, iist
         F_psi(1:sys%gr%mesh%np,idim) = zpsi(1:sys%gr%mesh%np)
      end do
      filename = 'floquet_td_spin_ik_'//trim(adjustl(str))//'_ist_'//trim(adjustl(str2))
      iunit = io_open(FLOQUET_DIR//filename, action='write')
      write(iunit,'(a)') '# time Sx Sy Sz'

      do it=1,hm%F%nT
        time = hm%F%Tcycle/hm%F%nT*it
        call floquet_td_state(hm%F,sys%gr%mesh,F_psi,time,psi_t)
        spin(1:3) = state_spin(sys%gr%mesh, psi_t) 
        write(iunit,'(e12.6,2x,e12.6,2x,e12.6,2x,e12.6)') time, spin(1:3)
      end do

      close(iunit)
    end do
  end do

  call restart_end(restart)

  SAFE_DEALLOCATE_A(psi_t)
  SAFE_DEALLOCATE_A(zpsi)
  SAFE_DEALLOCATE_A(F_psi)
  SAFE_DEALLOCATE_A(kpoints_input)
  SAFE_DEALLOCATE_A(states_input)
  end subroutine calc_floquet_td_spin

  !--------------------------------
  subroutine floquet_td_state(F,mesh,F_psi,time,psi_t)
  
  type(floquet_t) :: F
  type(mesh_t)    :: mesh
  CMPLX   :: F_psi(:,:)
  integer :: ik, ist
  FLOAT   :: time
  CMPLX   :: psi_t(:,:)
  
  integer :: idim, im, imm
  
  psi_t = M_ZERO
  do im= F%order(1),F%order(2)
    imm = im - F%order(1) + 1
    do idim=1,F%spindim
      psi_t(1:mesh%np,idim)  =  psi_t(1:mesh%np,idim) + exp(-M_zI*im*F%omega*time)*F_psi(1:mesh%np,(imm-1)*F%spindim+idim)
    end do
  end do

  ! normalize td-state
  psi_t(1:mesh%np,1:F%spindim)  = M_ONE/zmf_nrm2(mesh,F%spindim,psi_t)*psi_t(1:mesh%np,1:F%spindim) 

  end subroutine floquet_td_state
  
  !---------------------------------------
  subroutine floquet_photoelectron_spectrum(hm, sys, st, pomega, pol, spect, me)
    type(hamiltonian_t), intent(in) :: hm
    type(system_t), intent(in)      :: sys
    type(states_t), intent(in)      :: st
    FLOAT,          intent(in)      :: pomega     ! Probe field energy 
    FLOAT,          intent(in)      :: pol(:)     ! Probe field polarization vector
    FLOAT,          intent(out)     :: spect(:,:) ! the photoelectron spectrum
    FLOAT,          intent(out)     :: me(:,:)    ! the photoeletron matrix elements


    CMPLX, allocatable :: u_ma(:,:),  phase(:), tmp(:)
    FLOAT :: omega, dt , qq(1:3), kpt(1:3), xx(1:MAX_DIM)
    integer :: idx, im, it, ist, idim, nT, ik, ia, Fdim(2), imm, dim, pdim, ip, spindim

    type(mesh_t),   pointer :: mesh

    PUSH_SUB(floquet_photoelectron_spectrum)

    mesh  => sys%gr%der%mesh
    
    dt=hm%F%dt
    nT=hm%F%nT
    omega=hm%F%omega
    Fdim(:)=hm%F%order(:)
    spindim = hm%F%spindim 
    
    dim = mesh%sb%dim
    pdim = mesh%sb%periodic_dim
    
    SAFE_ALLOCATE(phase(1:mesh%np))
    SAFE_ALLOCATE(u_ma(1:mesh%np,hm%F%floquet_dim))
    SAFE_ALLOCATE(tmp(spindim))

    kpt(:) = M_ZERO
    qq(:)  = M_ZERO
    me(:,:) = M_ZERO
    spect(:,:) = M_ZERO
    do ik=st%d%kpt%start, st%d%kpt%end
      kpt(1:dim) = kpoints_get_point(mesh%sb%kpoints, ik) 

      do ia=st%st_start, st%st_end
        qq(dim)    = pomega - st%eigenval(ia,ik) - sum(kpt(1:pdim)**2)*M_HALF
        qq(1:pdim) = kpt(1:pdim)  
        
        do ip=1, mesh%np
          xx=mesh_x_global(mesh, ip) 
          phase(ip) = exp(-M_ZI*qq(dim)*xx(dim))
        end do
        
        call states_get_state(st, mesh, ia, ik, u_ma)
        
        tmp(:) = M_ZERO
        do idx=hm%F%flat_idx%start, hm%F%flat_idx%end
          it = hm%F%idx_map(idx,1)
          im = hm%F%idx_map(idx,2)
          imm = im - Fdim(1) + 1
          do idim=1,spindim
            tmp(idim)  =  tmp(idim) + exp(-M_zI*im*omega*it*dt)/nT * &
                          zmf_integrate(mesh, phase(1:mesh%np)*u_ma(1:mesh%np,(imm-1)*spindim+idim))
          end do
           
        end do 
        if(hm%F%is_parallel) call comm_allreduce(hm%F%mpi_grp%comm, tmp(:))   
        
        me(ia,ik) = sum(abs(tmp(:))**2) * sum((pol(1:dim)*qq(dim)))**2
        
        spect(ia,ik) =  me(ia,ik) * st%occ(ia,ik)
          
      end do
    end do
    
    call comm_allreduce(st%st_kpt_mpi_grp%comm, me)
    call comm_allreduce(st%st_kpt_mpi_grp%comm, spect)
    
    
    SAFE_DEALLOCATE_A(tmp)
    SAFE_DEALLOCATE_A(phase)
    SAFE_DEALLOCATE_A(u_ma)



    POP_SUB(floquet_photoelectron_spectrum)
    
  end subroutine floquet_photoelectron_spectrum    
  
  
  
  !--------------------------------------------
  subroutine calc_floquet_norms(mesh,kpoints,st,dressed_st,iter,floquet_dim)
    type(mesh_t), intent(in) :: mesh
    type(kpoints_t), intent(in) :: kpoints
    type(states_t), intent(in) :: st,dressed_st
    integer :: iter , floquet_dim

    integer :: maxiter, ik, in, im, ist, idim, ierr, nik, dim, nst, iunit
    CMPLX, allocatable :: temp_state1(:,:), temp_state2(:,:)
    FLOAT, allocatable :: norms(:,:,:)
    character(len=1024):: ik_name,iter_name, filename

    SAFE_ALLOCATE(temp_state1(1:mesh%np,st%d%dim))
    SAFE_ALLOCATE(temp_state2(1:mesh%np,floquet_dim*st%d%dim))
    SAFE_ALLOCATE(norms(1:kpoints%reduced%npoints,1:dressed_st%nst,floquet_dim))

    norms = M_ZERO
   

    do ik=st%d%kpt%start,st%d%kpt%end
      do in=1,floquet_dim
        do ist=st%st_start,st%st_end
          call states_get_state(dressed_st, mesh, (in-1)*st%nst+ist, ik,temp_state2)

          do im=1,floquet_dim
            do idim=1,st%d%dim
              temp_state1(1:mesh%np,idim) = temp_state2(1:mesh%np,(im-1)*st%d%dim+idim)
            end do
            norms(ik,(in-1)*st%nst+ist,im) = zmf_nrm2(mesh,st%d%dim,temp_state1)
          end do
        enddo
      enddo
    enddo
    
    call comm_allreduce(dressed_st%st_kpt_mpi_grp%comm, norms)

    if(mpi_world%rank==0) then
       write(iter_name,'(i4)') iter 
       do ik=kpoints%reduced%npoints-kpoints%nik_skip+1, kpoints%reduced%npoints
          write(ik_name,'(i4)') ik
          filename = FLOQUET_DIR//'/floquet_norms_ik_'//trim(adjustl(ik_name))//'_iter_'//trim(adjustl(iter_name))
          iunit = io_open(filename, action='write')
          
           do ist=1,dressed_st%nst
              do in=1,floquet_dim
                 write(iunit,'(i4,a1,e12.6,a1,i4,a1,e12.6)') ist, '', dressed_st%eigenval(ist,ik), ' ',in, ' ', norms(ik,ist,in)
             end do
             write(iunit,'(a1)') ' '
          end do

          call io_close(iunit)
       end do
    end if

    SAFE_DEALLOCATE_A(temp_state1)
    SAFE_DEALLOCATE_A(temp_state2)
    SAFE_DEALLOCATE_A(norms)
    
  end subroutine calc_floquet_norms
  

  end program floquet_observables

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
