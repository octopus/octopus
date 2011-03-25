
!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
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
subroutine PES_rc_init(pesrc, mesh, st, save_iter)
  type(PES_rc_t), intent(out) :: pesrc
  type(mesh_t),   intent(in)  :: mesh
  type(states_t), intent(in)  :: st
  integer,        intent(in)  :: save_iter

  type(block_t) :: blk
  integer  :: ip
  FLOAT :: xx(MAX_DIM)

  FLOAT :: dmin                                                                                                   
  integer :: rankmin
  
  PUSH_SUB(PES_rc_init)

  message(1) = 'Info: Calculating PES using rc technique.'
  call messages_info(1)

  !%Variable PES_rc_points
  !%Type block
  !%Section Time-Dependent::PES
  !%Description
  !% List of points at which to calculate the photoelectron spectrum by Suraud method.
  !% The exact syntax is:
  !%
  !% <tt>%PES_rc_points
  !% <br>&nbsp;&nbsp;x1 | y1 | z1
  !% <br>%
  !% </tt>
  !%End
  if (parse_block(datasets_check('PES_rc_points'), blk) < 0) call input_error('PES_rc_points')

  pesrc%npoints = parse_block_n(blk)

  ! setup filenames and read points
  SAFE_ALLOCATE(pesrc%filenames(1:pesrc%npoints))
  SAFE_ALLOCATE(pesrc%points   (1:pesrc%npoints))
  SAFE_ALLOCATE(pesrc%rankmin   (1:pesrc%npoints))

  do ip = 1, pesrc%npoints
    write(pesrc%filenames(ip), '(a,i2.2,a)') 'PES_rc.', ip, '.out'
    
    call parse_block_float(blk, ip - 1, 0, xx(1), units_inp%length)
    call parse_block_float(blk, ip - 1, 1, xx(2), units_inp%length)
    call parse_block_float(blk, ip - 1, 2, xx(3), units_inp%length)

    
    pesrc%points(ip) = mesh_nearest_point(mesh, xx, dmin, rankmin)
    pesrc%rankmin(ip)= rankmin

  end do

  call parse_block_end(blk)

  SAFE_ALLOCATE(pesrc%wf(1:pesrc%npoints, 1:st%d%dim, st%st_start:st%st_end, 1:st%d%nik, 1:save_iter))

  POP_SUB(PES_rc_init)
end subroutine PES_rc_init


! ---------------------------------------------------------
subroutine PES_rc_end(pesrc)
  type(PES_rc_t), intent(inout) :: pesrc

  PUSH_SUB(PES_rc_end)

  if(associated(pesrc%filenames)) then
    SAFE_DEALLOCATE_P(pesrc%filenames)
    SAFE_DEALLOCATE_P(pesrc%points)
    SAFE_DEALLOCATE_P(pesrc%wf)
    SAFE_DEALLOCATE_P(pesrc%rankmin)
  end if

  POP_SUB(PES_rc_end)
end subroutine PES_rc_end


! ---------------------------------------------------------
subroutine PES_rc_calc(pesrc, st,mesh, ii)
  type(PES_rc_t), intent(inout) :: pesrc
  type(states_t), intent(in)    :: st
  integer,        intent(in)    :: ii
  type(mesh_t),   intent(in) :: mesh

  integer :: ip, ik, ist, idim
  logical :: contains_ip
  complex(r8) :: wf
#if defined(HAVE_MPI)
  integer status(MPI_STATUS_SIZE)
#endif

  PUSH_SUB(PES_rc_calc)

  contains_ip = .true.

  do ip = 1, pesrc%npoints

#if defined(HAVE_MPI)
     if(mesh%mpi_grp%rank .eq. pesrc%rankmin(ip))then !needed if mesh%parallel_in_domains is true
        contains_ip = .true.
     else
        contains_ip = .false.
     end if
#endif
        
    do ik = st%d%kpt%start, st%d%kpt%end
      do ist = st%st_start, st%st_end
        do idim = 1, st%d%dim
           if(contains_ip) then
              pesrc%wf(ip, idim, ist, ik, ii) = st%occ(ist, ik) * &
                   st%zpsi(pesrc%points(ip), idim, ist, ik)

#if defined(HAVE_MPI)
              if(mesh%mpi_grp%rank .ne. 0) then
                 wf=pesrc%wf(ip, idim, ist, ik, ii)
                 call mpi_send(wf,1, MPI_DOUBLE_COMPLEX,0, 1, mesh%mpi_grp%comm, mpi_err)
              end if
#endif

           end if

#if defined(HAVE_MPI)
           if(mesh%mpi_grp%rank .eq. 0 .and. pesrc%rankmin(ip) .ne. 0) then
              call mpi_recv(wf,1, MPI_DOUBLE_COMPLEX,pesrc%rankmin(ip), 1, mesh%mpi_grp%comm,status, mpi_err)
              pesrc%wf(ip, idim, ist, ik, ii) = wf
              
           end if
#endif

        end do
      end do
    end do
  end do

  POP_SUB(PES_rc_calc)
end subroutine PES_rc_calc


! ---------------------------------------------------------
subroutine PES_rc_output(pesrc, st, iter, save_iter, dt)
  type(PES_rc_t), intent(in) :: pesrc
  type(states_t), intent(in) :: st
  integer,        intent(in) :: iter, save_iter
  FLOAT,          intent(in) :: dt

  type(unit_t) :: units

  integer :: ip, iunit, ii, jj, ik, ist, idim
  CMPLX :: vfu

  PUSH_SUB(PES_rc_output)

  do ip = 1, pesrc%npoints
    iunit = io_open('td.general/'//pesrc%filenames(ip), action='write', position='append')
    do ii = 1, save_iter
      jj = iter - save_iter + ii
      write(iunit, '(e17.10)', advance='no') units_from_atomic(units_inp%time, jj * dt)
      do ik = 1, st%d%nik
        do ist = st%st_start, st%st_end
          do idim = 1, st%d%dim
            vfu = pesrc%wf(ip, idim, ist, ik, ii)
            vfu = units_from_atomic(sqrt(units_out%length**(-3)), vfu)
            units=sqrt(units_out%length**(-3))
            write(iunit, '(1x,e18.10E3,1x,e18.10E3)', advance='no') &
                 real(vfu),  aimag(vfu) 
          end do
        end do
      end do
      write(iunit, '(1x)', advance='yes')
    end do
    call io_close(iunit)
  end do

  POP_SUB(PES_rc_output)
end subroutine PES_rc_output

! ---------------------------------------------------------
subroutine PES_rc_init_write(pesrc, mesh, st)
  type(PES_rc_t), intent(in) :: pesrc
  type(mesh_t),   intent(in)  :: mesh
  type(states_t), intent(in)  :: st

  integer  :: ip,ik,ist,idim,iunit

  FLOAT :: xx(MAX_DIM)



  PUSH_SUB(PES_rc_init_write)


  if(mpi_grp_is_root(mpi_world)) then

  
    do ip = 1, pesrc%npoints
      iunit = io_open('td.general/'//pesrc%filenames(ip), action='write')
      xx(:)=mesh%idx%Lxyz(pesrc%points(ip),:)
      write(iunit,'(a1)') '#'
      write(iunit, '(a7,f17.6,a1,f17.6,a1,f17.6,5a)') &
           '# R = (',units_from_atomic(units_inp%length, xx(1)*mesh%spacing(1)), &
           ' ,',units_from_atomic(units_inp%length, xx(2)*mesh%spacing(1)), &
           ' ,',units_from_atomic(units_inp%length, xx(3)*mesh%spacing(1)), &
             ' )  [', trim(units_abbrev(units_inp%length)), ']'

      write(iunit,'(a1)') '#'  
      write(iunit, '(a3,14x)', advance='no') '# t' 
      do ik = st%d%kpt%start, st%d%kpt%end
        do ist = st%st_start, st%st_end
          do idim = 1, st%d%dim
            write(iunit, '(3x,a8,i3,a7,i3,a8,i3,3x)', advance='no') &
                 "ik = ", ik, " ist = ", ist, " idim = ", idim
          end do
        end do
      end do
      write(iunit, '(1x)', advance='yes')
      
      call io_close(iunit)
    end do
    !!  end if

  endif




  POP_SUB(PES_rc_init_write)
end subroutine PES_rc_init_write


!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
