!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Electronic acceleration (to calculate harmonic spectrum...)
! It is calculated as:
!
! d2<x>/dt2 = d<p>/dt + i<[H,[V_nl,x]]> = 
!           = i<[V_l,p]> + i<[V_nl,p]> - E(t)N + i<[H,[V_nl,x]]>
!
! WARNING: This subroutine only works if ions are not allowed to move
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine td_calc_tacc(acc, t, reduce)
    implicit none
    FLOAT, intent(in)  :: t
    FLOAT, intent(out) :: acc(3)
    logical, intent(in), optional :: reduce

    FLOAT :: field(3), x(3), y(3), mesh_x(3), r, vl, dvl, d, charge
    CMPLX, allocatable :: hzpsi(:,:), hhzpsi(:,:), xzpsi(:,:,:), vnl_xzpsi(:,:)
    integer  :: j, k, is, i, ik, ist, idim, add_lm, l, m, ii, jj, ierr
    type(atom_type), pointer :: atm

    call push_sub('td_calc_tacc')

    ! The term i<[V_l,p]> + i<[V_nl,p]> may be considered as equal but opposite to the
    ! force exerted by the electrons on the ions. COMMENT: This has to be thought about.
    ! Maybe we are forgetting something....
    x = M_ZERO
    call zepot_forces(h%ep, sys)
    do i = 1, sys%natoms
      x = x - sys%atom(i)%f
    enddo
    acc = x
    
    ! Adds the laser contribution : i<[V_laser, p]>
    if(h%ep%no_lasers > 0) then
      call epot_laser_field(h%ep, t, field)
      acc(1:3) = acc(1:3) - sys%st%qtot*field(1:3)
    end if
    
    if(.not.sys%nlpp) then
      call pop_sub()
      return
    end if

    ! And now, i<[H,[V_nl,x]]>
    x = M_ZERO
    allocate(hzpsi(sys%m%np, sys%st%dim), hhzpsi(3, sys%m%np))

    do ik = 1, sys%st%nik
      do ist = sys%st%st_start, sys%st%st_end
        
        call zhpsi(h, sys%m, sys%st%zpsi(:, :, ist, ik), hzpsi(:,:), sys, ik)
        call epot_laser_field(h%ep, t, field)
        do k = 1, sys%m%np
          call mesh_xyz(sys%m, k, mesh_x)
          hzpsi(k,:) = hzpsi(k,:) + sum(mesh_x*field) * sys%st%zpsi(k,:,ist,ik)
        end do
        
        allocate(xzpsi(sys%m%np, sys%st%dim, 3), vnl_xzpsi(sys%m%np, sys%st%dim))
        xzpsi = M_z0
        do k = 1, sys%m%np
          do j = 1, conf%dim
            xzpsi(k, 1:sys%st%dim, j) = sys%m%Lxyz(j, k)*sys%m%h(j) * sys%st%zpsi(k, 1:sys%st%dim, ist, ik)
          end do
        end do
         
        do j = 1, conf%dim
          vnl_xzpsi = M_z0
          call zvnlpsi(h, sys%m, xzpsi(sys%m%np, 1:sys%st%dim, j), vnl_xzpsi(1:sys%m%np, 1:sys%st%dim), sys, ik)
               
          do idim = 1, sys%st%dim
            x(j) = x(j) - 2*sys%st%occ(ist, ik)*zmf_dotp(sys%m, R_CONJ(hzpsi(1:sys%m%np, idim)), &
                 vnl_xzpsi(1:sys%m%np, idim) )
          end do
        end do
           
        xzpsi = M_z0
        do k = 1, sys%m%np
          do j = 1, conf%dim
            xzpsi(k, 1:sys%st%dim, j) = sys%m%Lxyz(j, k)*sys%m%h(j) * hzpsi(k, 1:sys%st%dim)
          end do
        end do
        
        do j = 1, conf%dim
          vnl_xzpsi = M_z0
          call zvnlpsi(h, sys%m, xzpsi(sys%m%np, 1:sys%st%dim, j), vnl_xzpsi(1:sys%m%np, 1:sys%st%dim), sys, ik)
          do idim = 1, sys%st%dim
            x(j) = x(j) + 2*sys%st%occ(ist, ik)* &
                  zmf_dotp(sys%m, R_CONJ(sys%st%zpsi(1:sys%m%np, idim, ist, ik)), &
                  vnl_xzpsi(1:sys%m%np, idim) )
          end do
        end do
        deallocate(xzpsi, vnl_xzpsi)
       
     end do
   end do
   deallocate(hzpsi, hhzpsi)

#if defined(HAVE_MPI) && defined(MPI_TD)
   if(present(reduce)) then
     if(reduce) then
       call MPI_ALLREDUCE(x(1), y(1), conf%dim, &
            MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
       x = y
     end if
   end if
#endif
   acc = acc + x
   
   call pop_sub()
 end subroutine td_calc_tacc
