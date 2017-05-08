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
  use messages_oct_m
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
  
!   type(simul_box_t)    :: sb
!   type(grid_t)         :: gr
  type(restart_t)      :: restart
  
  type(system_t)      :: sys
    
  type(hamiltonian_t) :: hm
  
  character(len=512)   :: filename, str

    
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
  call calc_mode_par_set_parallelization(P_STRATEGY_OTHER, default = .true.)

  call fft_all_init()
  call unit_system_init()
  
  call system_init(sys)


  call floquet_init(sys,hm%F,sys%st%d%dim)
  gs_st => sys%st

  call states_init(dressed_st, sys%gr, sys%geo,floquet_dim=hm%F%floquet_dim)
  call kpoints_distribute(dressed_st%d,sys%mc)
  call states_distribute_nodes(dressed_st,sys%mc)
  call states_exec_init(dressed_st, sys%mc)
  call states_allocate_wfns(dressed_st,sys%gr%der%mesh, wfs_type = TYPE_CMPLX)


  
  
  call restart_module_init()
  
  call floquet_restart_dressed_st(hm, sys, dressed_st, ierr)
  
  
  call messages_write('Read Floquet restart files.')
  call messages_info()

  !%Variable FloquetObservableCalc
  !%Type flag
  !%Default none
  !%Section Utilities::oct-floquet_observables
  !%Description
  !% Specifies what observables to calculate form the Floquet solution.
  !% Example: <tt>f_norms + f_spin</tt>
  !%Option f_norms bit(1)
  !% The norms of the Floquet states with each tensorial subspace.
  !%Option f_arpes bit(2)
  !% Calculate ARPES matrix elements for Floquet states.
  !%Option f_spin bit(3)
  !% Calculate the spin polarization of each state. 
  !%End
  call parse_variable('FloquetObservableCalc', out_what, out_what)
  
!   call getopt_floquet_observables(uEstep, uEspan,&
!                                      uThstep, uThspan, uPhstep, &
!                                      uPhspan, pol, center, pvec, integrate)
!
!   call getopt_end()
                                       

  if(iand(out_what, OPTION__FLOQUETOBSERVABLECALC__F_NORMS) /= 0) then
    call messages_write('Calculate norms of Floquet subspaces.')
    call messages_info()
    call floquet_calc_norms(sys%gr%der%mesh,sys%gr%sb%kpoints,gs_st,dressed_st, hm%F%iter,hm%F%floquet_dim)
  end if

  if(iand(out_what, OPTION__FLOQUETOBSERVABLECALC__F_ARPES) /= 0) then
    call messages_write('Calculate Floqeut ARPES.')
    call messages_info()

    call calc_floquet_arpes()
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
  


end program floquet_observables

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
