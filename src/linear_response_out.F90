!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch, X. Andrade.
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
!! -*- coding: utf-8 mode: f90 -*-
!! $Id$

! ---------------------------------------------------------
subroutine X(lr_output) (st, gr, lr, dir, tag, outp)
  type(states_t),   intent(inout) :: st
  type(grid_t),     intent(inout) :: gr
  type(lr_t),       intent(inout) :: lr
  character(len=*), intent(in)    :: dir
  integer,          intent(in)    :: tag
  type(output_t),   intent(in)    :: outp

  integer :: ik, ist, idim, ierr, is, i
  character(len=80) :: fname
  FLOAT :: u
  FLOAT, allocatable :: dtmp(:)
  R_TYPE, allocatable :: tmp(:)
  

  call push_sub('lr_response_out.lr_output')

  u = M_ONE/units_out%length%factor**NDIM

  if(iand(outp%what, output_density).ne.0) then
    do is = 1, st%d%nspin
      write(fname, '(a,i1,a,i1)') 'lr_density-', is, '-', tag
      call X(output_function)(outp%how, dir, fname, gr%m, gr%sb, lr%X(dl_rho)(:, is), u, ierr)
    end do
  end if

  if(iand(outp%what, output_pol_density).ne.0) then
    ALLOCATE(tmp(1:NP),NP)
    do is = 1, st%d%nspin
      do i=1,NDIM
        tmp(1:NP)=gr%m%x(1:NP,i)*lr%X(dl_rho)(:, is)
        write(fname, '(a,i1,a,i1,a,i1)') 'alpha_density-', is, '-', tag, '-', i
        call X(output_function)(outp%how, dir, fname, gr%m, gr%sb, tmp, u, ierr)
      end do
    end do
    deallocate(tmp)
  end if

  if( (iand(outp%what, output_current).ne.0) .and. (st%d%wfs_type == M_CMPLX) )then
    do is = 1, st%d%nspin
      do idim = 1, NDIM
        write(fname, '(a,i1,a,i1,a,a)') 'lr_current-', is, '-', tag, '-',  index2axis(idim)
        call zoutput_function(outp%how, dir, fname, gr%m, gr%sb, lr%dl_j(:, idim, is), u, ierr)
      end do
    end do
  end if


  if(iand(outp%what, output_wfs).ne.0) then
    do ist = st%st_start, st%st_end
      if(loct_isinstringlist(ist, outp%wfs_list)) then
        do ik = 1, st%d%nik
          do idim = 1, st%d%dim
            write(fname, '(a,i3.3,a,i3.3,a,i1,a,i1)') 'lr_wf-', ik, '-', ist, '-', idim,'-', tag
            call X(output_function) (outp%how, dir, fname, gr%m, gr%sb, &
              lr%X(dl_psi) (1:, idim, ist, ik), sqrt(u), ierr)
          end do
        end do
      end if
    end do
  end if

  if(iand(outp%what, output_wfs_sqmod).ne.0) then
    ALLOCATE(dtmp(NP_PART), NP_PART)
    do ist = 1, st%nst
      if(loct_isinstringlist(ist, outp%wfs_list)) then
        do ik = 1, st%d%nik
          do idim = 1, st%d%dim
            write(fname, '(a,i3.3,a,i3.3,a,i1)') 'sqm_lr_wf-', ik, '-', ist, '-', idim
            dtmp = abs(lr%X(dl_psi) (:, idim, ist, ik))**2
            call doutput_function (outp%how, dir, fname, gr%m, gr%sb, dtmp, u, ierr)
          end do
        end do
      end if
    end do
    deallocate(dtmp)
  end if

  if(NDIM==3) then
    if(iand(outp%what, output_elf).ne.0)    call lr_elf('lr_D','lr_elf')
  end if

  call pop_sub()
contains

  ! ---------------------------------------------------------
  subroutine lr_elf(filename1, filename2)
    character(len=*), intent(in) :: filename1
    character(len=*), intent(in) :: filename2
    
    integer :: is, ierr
    u=M_ONE

    do is = 1, st%d%nspin
      write(fname, '(a,a,i1,a,i1)') trim(filename1), '-', is, '-', tag 
      call X(output_function)(outp%how, dir, trim(fname), gr%m, gr%sb, lr%X(dl_de)(1:NP,is), u, ierr)
    end do

    do is = 1, st%d%nspin
      write(fname, '(a,a,i1,a,i1)') trim(filename2), '-', is, '-', tag 
      call X(output_function)(outp%how, dir, trim(fname), gr%m, gr%sb, lr%X(dl_elf)(1:NP,is), u, ierr)
    end do

  end subroutine lr_elf

end subroutine X(lr_output)

