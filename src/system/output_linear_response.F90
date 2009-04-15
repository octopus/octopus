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
!! $Id$

! ---------------------------------------------------------
subroutine X(h_sys_output_lr) (st, gr, lr, dir, tag, isigma, outp, geo)
  type(states_t),       intent(inout) :: st
  type(grid_t),         intent(inout) :: gr
  type(lr_t),           intent(inout) :: lr
  character(len=*),     intent(in)    :: dir
  integer,              intent(in)    :: tag, isigma
  type(h_sys_output_t), intent(in)    :: outp
  type(geometry_t),     intent(in)    :: geo

  integer :: ik, ist, idim, ierr, is, i
  character(len=80) :: fname
  FLOAT :: u
  FLOAT, allocatable :: dtmp(:)
  R_TYPE, allocatable :: tmp(:)
  

  call push_sub('output_linear_response.Xh_sys_output_lr')

  u = M_ONE/units_out%length%factor**gr%mesh%sb%dim

  if(isigma == 1) then ! the density, current, etc. are only defined for the + frequency

    if(iand(outp%what, output_density).ne.0) then
      do is = 1, st%d%nspin
        write(fname, '(a,i1,a,i1)') 'lr_density-', is, '-', tag
        call X(output_function)(outp%how, dir, fname, gr%mesh, gr%sb, lr%X(dl_rho)(:, is), u, ierr, geo = geo)
      end do
    end if

    if(iand(outp%what, output_pol_density).ne.0) then
      ALLOCATE(tmp(1:gr%mesh%np),gr%mesh%np)
      do is = 1, st%d%nspin
        do i=1,gr%mesh%sb%dim
          tmp(1:gr%mesh%np)=gr%mesh%x(1:gr%mesh%np,i)*lr%X(dl_rho)(:, is)
          write(fname, '(a,i1,a,i1,a,i1)') 'alpha_density-', is, '-', tag, '-', i
          call X(output_function)(outp%how, dir, fname, gr%mesh, gr%sb, tmp, u, ierr, geo = geo)
        end do
      end do
      SAFE_DEALLOCATE_A(tmp)
    end if

    if( (iand(outp%what, output_current).ne.0) .and. (st%wfs_type == M_CMPLX) )then
      do is = 1, st%d%nspin
        do idim = 1, gr%mesh%sb%dim
          write(fname, '(a,i1,a,i1,a,a)') 'lr_current-', is, '-', tag, '-',  index2axis(idim)
          call zoutput_function(outp%how, dir, fname, gr%mesh, gr%sb, lr%dl_j(:, idim, is), u, ierr, geo = geo)
        end do
      end do
    end if

    if(gr%mesh%sb%dim==3) then
      if(iand(outp%what, output_elf).ne.0) call lr_elf('lr_D','lr_elf')
    end if

  end if ! isigma == 1


  if(iand(outp%what, output_wfs).ne.0) then
    do ist = st%st_start, st%st_end
      if(loct_isinstringlist(ist, outp%wfs_list)) then
        do ik = st%d%kpt%start, st%d%kpt%end
          do idim = 1, st%d%dim
            write(fname, '(a,i3.3,a,i3.3,a,i1,a,i1,a,i1)') &
              'lr_wf-', ik, '-', ist, '-', idim, '-', tag, '-', isigma
            call X(output_function) (outp%how, dir, fname, gr%mesh, gr%sb, &
              lr%X(dl_psi) (1:, idim, ist, ik), sqrt(u), ierr, geo = geo)
          end do
        end do
      end if
    end do
  end if

  if(iand(outp%what, output_wfs_sqmod).ne.0) then
    ALLOCATE(dtmp(gr%mesh%np_part), gr%mesh%np_part)
    do ist = 1, st%nst
      if(loct_isinstringlist(ist, outp%wfs_list)) then
         do ik = st%d%kpt%start, st%d%kpt%end
          do idim = 1, st%d%dim
            write(fname, '(a,i3.3,a,i3.3,a,i1,a,i1,a,i1)') &
              'sqm_lr_wf-', ik, '-', ist, '-', idim, '-', isigma

            dtmp = abs(lr%X(dl_psi) (:, idim, ist, ik))**2
            call doutput_function (outp%how, dir, fname, gr%mesh, gr%sb, dtmp, u, ierr, geo = geo)
          end do
        end do
      end if
    end do
    SAFE_DEALLOCATE_A(dtmp)
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
      call X(output_function)(outp%how, dir, trim(fname), gr%mesh, gr%sb, lr%X(dl_de)(1:gr%mesh%np,is), u, ierr, geo = geo)
    end do

    do is = 1, st%d%nspin
      write(fname, '(a,a,i1,a,i1)') trim(filename2), '-', is, '-', tag 
      call X(output_function)(outp%how, dir, trim(fname), gr%mesh, gr%sb, lr%X(dl_elf)(1:gr%mesh%np,is), u, ierr, geo = geo)
    end do

  end subroutine lr_elf

end subroutine X(h_sys_output_lr)


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
