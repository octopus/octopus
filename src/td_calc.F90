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

  subroutine td_calc_projection(p)
    complex(r4), intent(out) :: p(u_st%nst, sys%st%st_start:sys%st%st_end, sys%st%nik)

    integer :: uist, uik, ist, ik

    do ik = 1, sys%st%nik
      do ist = sys%st%st_start, sys%st%st_end
        do uist = 1, u_st%nst
          p(uist, ist, ik) = cmplx(sum(sys%st%zpsi(1:sys%m%np,:, ist, ik)* &
               u_st%R_FUNC(psi) (1:sys%m%np,:, uist, ik)), kind=r4)*sys%m%vol_pp
        end do
      end do
    end do
  end subroutine td_calc_projection

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Electronic acceleration (to calculate harmonic spectrum...)
! It is calculated as:
!
! d2<x>/dt2 = d<p>/dt + i<[H,[V_nl,x]]> = 
!           = i<[V_l,p]> + i<[V_nl,p]> - E(t)N + i<[H,[V_nl,x]]>
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine td_calc_tacc(acc, t, reduce)
    implicit none
    real(r8), intent(in)  :: t
    real(r8), intent(out) :: acc(3)
    logical, intent(in), optional :: reduce

    real(r8) :: field(3), x(3), y(3), mesh_x(3), r, vl, dvl, d, charge
    complex(r8), allocatable :: hzpsi(:,:), hhzpsi(:,:), xzpsi(:,:,:), vnl_xzpsi(:,:)
    integer  :: j, k, is, i, ik, ist, idim, add_lm, l, m, ii, jj, ierr
    type(atom_type), pointer :: atm

    sub_name = 'td_calc_tacc'; call push_sub()

    ! This to make sure ions do not move...
    if(td%move_ions > 0) then
      message(1) = 'Error. If harmonic spectrum is to be calculated, moves should not move'
      message(2) = '(In present version)'
      call write_fatal(2)
    end if

    ! The term i<[V_l,p]> + i<[V_nl,p]> may be considered as equal but opposite to the
    ! force exerted by the electrons on the ions. COMMENT: This has to be thought about.
    ! Maybe we are forgetting something....
    x = 0.0_r8
    call zforces(h, sys)
    do i = 1, sys%natoms
       x = x - sys%atom(i)%f
    enddo
    acc = x

    ! Adds the laser contribution : i<[V_laser, p]>
    if(td%no_lasers > 0) then
      call laser_field(td%no_lasers, td%lasers, t, field)
      acc(1:3) = acc(1:3) - sys%st%qtot*field(1:3)
    end if

    ! And now, i<[H,[V_nl,x]]>
    x = 0.0_r8
    allocate(xzpsi(0:sys%m%np, sys%st%dim, 3), vnl_xzpsi(sys%m%np, sys%st%dim), &
             hzpsi(sys%m%np, sys%st%dim), hhzpsi(3, sys%m%np) )
    do ik = 1, sys%st%nik
       do ist = sys%st%st_start, sys%st%st_end

            call zhpsi(h, sys, ik, sys%st%zpsi(:, :, ist, ik), hzpsi(:,:))
            call laser_field(td%no_lasers, td%lasers, t, field)
            do k = 1, sys%m%np
              call mesh_xyz(sys%m, k, mesh_x)
              hzpsi(k,:) = hzpsi(k,:) + sum(mesh_x*field) * sys%st%zpsi(k,:,ist,ik)
            end do

            xzpsi = (0.0_r8, 0.0_r8)
            do k = 1, sys%m%np
               xzpsi(k, 1:sys%st%dim, 1) = sys%m%lx(k)*sys%m%h(1) * sys%st%zpsi(k, 1:sys%st%dim, ist, ik)
               xzpsi(k, 1:sys%st%dim, 2) = sys%m%ly(k)*sys%m%h(2) * sys%st%zpsi(k, 1:sys%st%dim, ist, ik)
               xzpsi(k, 1:sys%st%dim, 3) = sys%m%lz(k)*sys%m%h(3) * sys%st%zpsi(k, 1:sys%st%dim, ist, ik)
            enddo

            do j = 1, 3
               vnl_xzpsi = (0.0_r8, 0.0_r8)
               call zvnlpsi( sys, xzpsi(0:sys%m%np, 1:sys%st%dim, j), vnl_xzpsi(1:sys%m%np, 1:sys%st%dim))
               do idim = 1, sys%st%dim
                  x(j) = x(j) - 2*sys%st%occ(ist, ik)*zmesh_dotp(sys%m, R_CONJ(hzpsi(1:sys%m%np, idim)), &
                                                                   vnl_xzpsi(1:sys%m%np, idim) )
               enddo
            enddo

            xzpsi = (0.0_r8, 0.0_r8)
            do k = 1, sys%m%np
               xzpsi(k, 1:sys%st%dim, 1) = sys%m%lx(k)*sys%m%h(1) * hzpsi(k, 1:sys%st%dim)
               xzpsi(k, 1:sys%st%dim, 2) = sys%m%ly(k)*sys%m%h(2) * hzpsi(k, 1:sys%st%dim)
               xzpsi(k, 1:sys%st%dim, 3) = sys%m%lz(k)*sys%m%h(3) * hzpsi(k, 1:sys%st%dim)
            enddo

            do j = 1, 3
               vnl_xzpsi = (0.0_r8, 0.0_r8)
               call zvnlpsi( sys, xzpsi(0:sys%m%np, 1:sys%st%dim, j), vnl_xzpsi(1:sys%m%np, 1:sys%st%dim))
               do idim = 1, sys%st%dim
                  x(j) = x(j) + 2*sys%st%occ(ist, ik)* zmesh_dotp(sys%m, R_CONJ(sys%st%zpsi(1:sys%m%np, idim, ist, ik)), &
                                            vnl_xzpsi(1:sys%m%np, idim) )
               enddo
            enddo

        enddo
    enddo
    deallocate(xzpsi, hzpsi, hhzpsi, vnl_xzpsi)

#if defined(HAVE_MPI) && defined(MPI_TD)
    if(present(reduce) .and. reduce) then
       call MPI_ALLREDUCE(x(1), y(1), 3, &
                          MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
       x = y
    end if
#endif
    acc = acc + x

    call pop_sub(); return 
  end subroutine td_calc_tacc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! APPENDIX.
!
! Next commented subroutines are for developing/testing purposes. They have
! the same interface as "td_calc_tacc". They calculate either the velocity
! or the acceleration in various ways. Should be commented in a not-developing
! version of the code.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  ! d<x>/dt = <p> + i<[V_nl, x]>
!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  subroutine td_calc_vel1(acc, t, reduce)
!!$    implicit none
!!$    real(r8), intent(in)  :: t
!!$    real(r8), intent(out) :: acc(3)
!!$    logical, intent(in), optional :: reduce
!!$
!!$    real(r8) :: x(3), y(3), d
!!$    complex(r8) :: p, z(3)
!!$    integer  :: j, k, is, i, ik, ist, idim, add_lm, l, m, ii, jj, ierr
!!$    complex(r8), allocatable :: hzpsi(:,:), hhzpsi(:,:), xzpsi(:,:,:), vnl_xzpsi(:,:)
!!$
!!$    sub_name = 'td_calc_vel1'; call push_sub()
!!$
!!$    x = 0.0_r8
!!$    ! This calculates <p>
!!$    allocate(hzpsi(3, sys%m%np))
!!$    do ik = 1, sys%st%nik
!!$       do ist = sys%st%st_start, sys%st%st_end
!!$          do idim = 1, sys%st%dim
!!$             call zmesh_derivatives (sys%m, sys%st%zpsi(:, idim, ist, ik), grad = hzpsi, &
!!$                                  alpha = -M_zI )
!!$             do j = 1, 3
!!$                x(j) = x(j) + &
!!$                 sys%st%occ(ist,ik)*zmesh_dotp(sys%m, R_CONJ(sys%st%zpsi(1:sys%m%np, idim,ist,ik)), &
!!$                                               hzpsi(j,1:sys%m%np) )
!!$             enddo
!!$          enddo
!!$       enddo
!!$    enddo
!!$    deallocate(hzpsi)
!!$
!!$    ! And this calculates i<[V_nl, x]>
!!$    allocate(xzpsi(0:sys%m%np, sys%st%dim, 3), vnl_xzpsi(sys%m%np, sys%st%dim))
!!$    do ik = 1, sys%st%nik
!!$       do ist = sys%st%st_start, sys%st%st_end
!!$
!!$            xzpsi(0, :, :) = (0.0_r8, 0.0_r8)
!!$            do k = 1, sys%m%np
!!$               xzpsi(k, 1:sys%st%dim, 1) = sys%m%lx(k)*sys%m%h(1) * sys%st%zpsi(k, 1:sys%st%dim, ist, ik)
!!$               xzpsi(k, 1:sys%st%dim, 2) = sys%m%ly(k)*sys%m%h(2) * sys%st%zpsi(k, 1:sys%st%dim, ist, ik)
!!$               xzpsi(k, 1:sys%st%dim, 3) = sys%m%lz(k)*sys%m%h(3) * sys%st%zpsi(k, 1:sys%st%dim, ist, ik)
!!$            enddo
!!$
!!$            do j = 1, 3
!!$               vnl_xzpsi = (0.0_r8, 0.0_r8)
!!$               call zvnlpsi( sys, xzpsi(0:sys%m%np, 1:sys%st%dim, j), vnl_xzpsi(1:sys%m%np, 1:sys%st%dim))
!!$               p = (0.0_r8, 0.0_r8)
!!$               do idim = 1, sys%st%dim
!!$                  p = p + sys%st%occ(ist, ik)*zmesh_dotp(sys%m, R_CONJ(sys%st%zpsi(1:sys%m%np, idim, ist, ik)), &
!!$                                              vnl_xzpsi(1:sys%m%np, idim) )
!!$               enddo
!!$               x(j) = x(j) - 2.0_r8*aimag(p)
!!$            enddo
!!$
!!$        enddo
!!$    enddo
!!$
!!$#if defined(HAVE_MPI) && defined(MPI_TD)
!!$    if(present(reduce) .and. reduce) then
!!$       call MPI_ALLREDUCE(x(1), y(1), 3, &
!!$                          MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
!!$       x = y
!!$    end if
!!$#endif
!!$    acc = x
!!$    deallocate(xzpsi, vnl_xzpsi)
!!$
!!$    call pop_sub(); return
!!$  end subroutine td_calc_vel1


!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  ! d<x>/dt = i<[H, x]>
!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  subroutine td_calc_vel2(acc, t, reduce)
!!$    implicit none
!!$    real(r8), intent(in)  :: t
!!$    real(r8), intent(out) :: acc(3)
!!$    logical, intent(in), optional :: reduce
!!$
!!$    integer :: ik, ist, k, j, idim, ierr
!!$    real(r8) :: x(3), y(3), field(3), mesh_x(3)
!!$    complex(r8) :: p
!!$    complex(r8), allocatable :: xzpsi(:,:,:), hzpsi(:,:)
!!$
!!$    sub_name = 'td_calc_tacc'; call push_sub()
!!$
!!$    x = 0.0_r8
!!$    allocate(xzpsi(0:sys%m%np, sys%st%dim, 3), hzpsi(sys%m%np, sys%st%dim))
!!$    do ik = 1, sys%st%nik
!!$       do ist = sys%st%st_start, sys%st%st_end
!!$
!!$            hzpsi = (0.0_r8, 0.0_r8)
!!$            call zhpsi(h, sys, ik, sys%st%zpsi(:, :, ist, ik), hzpsi(:,:))
!!$            call laser_field(td%no_lasers, td%lasers, t, field)
!!$            do k = 1, sys%m%np
!!$              call mesh_xyz(sys%m, k, mesh_x)
!!$              hzpsi(k,:) = hzpsi(k,:) + sum(mesh_x*field) * sys%st%zpsi(k,:,ist,ik)
!!$            end do
!!$            hzpsi = R_CONJ( hzpsi )
!!$
!!$            xzpsi(0, :, :) = (0.0_r8, 0.0_r8)
!!$            do k = 1, sys%m%np
!!$               xzpsi(k, 1:sys%st%dim, 1) = sys%m%lx(k)*sys%m%h(1) * sys%st%zpsi(k, 1:sys%st%dim, ist, ik)
!!$               xzpsi(k, 1:sys%st%dim, 2) = sys%m%ly(k)*sys%m%h(2) * sys%st%zpsi(k, 1:sys%st%dim, ist, ik)
!!$               xzpsi(k, 1:sys%st%dim, 3) = sys%m%lz(k)*sys%m%h(3) * sys%st%zpsi(k, 1:sys%st%dim, ist, ik)
!!$            enddo
!!$
!!$            do j = 1, 3
!!$               p = (0.0_r8, 0.0_r8)
!!$               do idim = 1, sys%st%dim
!!$                  p = p + sys%st%occ(ist, ik) * zmesh_dotp(sys%m, hzpsi(:, idim), xzpsi(1:, idim, j))
!!$               enddo
!!$               x(j) = x(j) - 2.0_r8*aimag(p)
!!$            enddo
!!$
!!$        enddo
!!$    enddo
!!$    deallocate(xzpsi, hzpsi)
!!$
!!$#if defined(HAVE_MPI) && defined(MPI_TD)
!!$    if(present(reduce) .and. reduce) then
!!$       call MPI_ALLREDUCE(x(1), y(1), 3, &
!!$                          MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
!!$       x = y
!!$    end if
!!$#endif
!!$    acc = x
!!$
!!$    call pop_sub(); return
!!$  end subroutine td_calc_vel2


!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  ! d2<x>/dt2 = -<[H,[H, x]]>
!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  subroutine td_calc_tacc1(acc, t, reduce)
!!$    implicit none
!!$    real(r8), intent(in)  :: t
!!$    real(r8), intent(out) :: acc(3)
!!$    logical, intent(in), optional :: reduce
!!$
!!$    real(r8) :: field(3), x(3), y(3), mesh_x(3)
!!$    real(r8), save :: prev(1:3) = 0.0_r8
!!$    integer  :: j, k, is, i, ik, ist, idim, ierr
!!$
!!$    complex(r8), allocatable :: hzpsi(:,:), hhzpsi(:,:), xzpsi(:,:,:)
!!$
!!$    sub_name = 'td_calc_tacc'; call push_sub()
!!$
!!$    x = 0.0_r8
!!$    allocate(xzpsi(0:sys%m%np, sys%st%dim, 3), hzpsi(0:sys%m%np, sys%st%dim), hhzpsi(0:sys%m%np, sys%st%dim))
!!$    do ik = 1, sys%st%nik
!!$       do ist = sys%st%st_start, sys%st%st_end
!!$
!!$            hzpsi(0:sys%m%np, 1:sys%st%dim) = (0.0_r8, 0.0_r8)
!!$            hhzpsi(0:sys%m%np, 1:sys%st%dim) = (0.0_r8, 0.0_r8)
!!$            call zhpsi(h, sys, ik, sys%st%zpsi(:, :, ist, ik), hzpsi(1:,:))
!!$            call laser_field(td%no_lasers, td%lasers, t, field)
!!$            do k = 1, sys%m%np
!!$              call mesh_xyz(sys%m, k, mesh_x)
!!$              hzpsi(k,:) = hzpsi(k,:) + sum(mesh_x*field) * sys%st%zpsi(k,:,ist,ik)
!!$            end do
!!$            call zhpsi(h, sys, ik, hzpsi(:, :), hhzpsi(1:,:))
!!$            do k = 1, sys%m%np
!!$              call mesh_xyz(sys%m, k, mesh_x)
!!$              hhzpsi(k,:) = hhzpsi(k,:) + sum(mesh_x*field) * hzpsi(k,:)
!!$            end do
!!$            hhzpsi = R_CONJ(hhzpsi) ! Warning: this has been complex conjugated here....
!!$
!!$            xzpsi = (0.0_r8, 0.0_r8)
!!$            do k = 1, sys%m%np
!!$               xzpsi(k, 1:sys%st%dim, 1) = sys%m%lx(k)*sys%m%h(1) * sys%st%zpsi(k, 1:sys%st%dim, ist, ik)
!!$               xzpsi(k, 1:sys%st%dim, 2) = sys%m%ly(k)*sys%m%h(2) * sys%st%zpsi(k, 1:sys%st%dim, ist, ik)
!!$               xzpsi(k, 1:sys%st%dim, 3) = sys%m%lz(k)*sys%m%h(3) * sys%st%zpsi(k, 1:sys%st%dim, ist, ik)
!!$            enddo
!!$
!!$            do j = 1, 3
!!$               do idim = 1, sys%st%dim
!!$                  x(j) = x(j) - 2.0_r8*sys%st%occ(ist, ik)*zmesh_dotp(sys%m, hhzpsi(1:, idim), xzpsi(1:, idim, j))
!!$               enddo
!!$            enddo
!!$
!!$            xzpsi = (0.0_r8, 0.0_r8)
!!$            do k = 1, sys%m%np
!!$               xzpsi(k, 1:sys%st%dim, 1) = sys%m%lx(k)*sys%m%h(1) * hzpsi(k, 1:sys%st%dim)
!!$               xzpsi(k, 1:sys%st%dim, 2) = sys%m%ly(k)*sys%m%h(2) * hzpsi(k, 1:sys%st%dim)
!!$               xzpsi(k, 1:sys%st%dim, 3) = sys%m%lz(k)*sys%m%h(3) * hzpsi(k, 1:sys%st%dim)
!!$            enddo
!!$            hzpsi = R_CONJ(hzpsi)  ! And also this is complex conjugated here...
!!$
!!$            do j = 1, 3
!!$               do idim = 1, sys%st%dim
!!$                  x(j) = x(j) + 2.0_r8*sys%st%occ(ist, ik)*zmesh_dotp(sys%m, hzpsi(1:, idim), xzpsi(1:, idim, j))
!!$               enddo
!!$            enddo
!!$
!!$        enddo
!!$    enddo
!!$
!!$#if defined(HAVE_MPI) && defined(MPI_TD)
!!$    if(present(reduce) .and. reduce) then
!!$       call MPI_ALLREDUCE(x(1), y(1), 3, &
!!$                          MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
!!$       x = y
!!$    end if
!!$#endif
!!$    acc = x
!!$    deallocate(xzpsi, hzpsi, hhzpsi)
!!$
!!$    call pop_sub(); return
!!$  end subroutine td_calc_tacc1


!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  ! d2<x>/dt2 = i<[H, p]> - <[H,[H,V_nl]]>
!!$  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  subroutine td_calc_tacc2(acc, t, reduce)
!!$    implicit none
!!$    real(r8), intent(in)  :: t
!!$    real(r8), intent(out) :: acc(3)
!!$    logical, intent(in), optional :: reduce
!!$
!!$    real(r8) :: field(3), x(3), y(3), mesh_x(3), r, d
!!$    integer  :: j, k, is, i, ik, ist, idim, ierr
!!$    complex(r8), allocatable :: hzpsi(:,:), xzpsi(:,:,:), vnl_xzpsi(:,:), hhzpsi(:,:)
!!$
!!$    sub_name = 'td_calc_tacc'; call push_sub()
!!$
!!$    x = 0.0_r8
!!$    allocate(xzpsi(0:sys%m%np, sys%st%dim, 3), vnl_xzpsi(sys%m%np, sys%st%dim), &
!!$             hzpsi(sys%m%np, sys%st%dim), hhzpsi(3, sys%m%np) )
!!$
!!$    do ik = 1, sys%st%nik
!!$       do ist = sys%st%st_start, sys%st%st_end
!!$
!!$            call zhpsi(h, sys, ik, sys%st%zpsi(:, :, ist, ik), hzpsi(:,:))
!!$            call laser_field(td%no_lasers, td%lasers, t, field)
!!$            do k = 1, sys%m%np
!!$              call mesh_xyz(sys%m, k, mesh_x)
!!$              hzpsi(k,:) = hzpsi(k,:) + sum(mesh_x*field) * sys%st%zpsi(k,:,ist,ik)
!!$            end do
!!$
!!$            do j = 1, 3
!!$               do idim = 1, sys%st%dim
!!$                  call zmesh_derivatives (sys%m, sys%st%zpsi(:, idim, ist, ik), grad = hhzpsi, alpha = -M_zI)
!!$                  x(j) = x(j) - 2*aimag( sys%st%occ(ist, ik)*zmesh_dotp(sys%m, R_CONJ(hzpsi(1:,idim)), hhzpsi(j, 1:)) )
!!$               enddo
!!$            enddo
!!$
!!$            xzpsi = (0.0_r8, 0.0_r8)
!!$            do k = 1, sys%m%np
!!$               xzpsi(k, 1:sys%st%dim, 1) = sys%m%lx(k)*sys%m%h(1) * sys%st%zpsi(k, 1:sys%st%dim, ist, ik)
!!$               xzpsi(k, 1:sys%st%dim, 2) = sys%m%ly(k)*sys%m%h(2) * sys%st%zpsi(k, 1:sys%st%dim, ist, ik)
!!$               xzpsi(k, 1:sys%st%dim, 3) = sys%m%lz(k)*sys%m%h(3) * sys%st%zpsi(k, 1:sys%st%dim, ist, ik)
!!$            enddo
!!$
!!$            do j = 1, 3
!!$               vnl_xzpsi = (0.0_r8, 0.0_r8)
!!$               call zvnlpsi( sys, xzpsi(0:sys%m%np, 1:sys%st%dim, j), vnl_xzpsi(1:sys%m%np, 1:sys%st%dim))
!!$               do idim = 1, sys%st%dim
!!$                  x(j) = x(j) - 2*sys%st%occ(ist, ik)*zmesh_dotp(sys%m, R_CONJ(hzpsi(1:sys%m%np, idim)), &
!!$                                                                   vnl_xzpsi(1:sys%m%np, idim) )
!!$               enddo
!!$            enddo
!!$
!!$            xzpsi = (0.0_r8, 0.0_r8)
!!$            do k = 1, sys%m%np
!!$               xzpsi(k, 1:sys%st%dim, 1) = sys%m%lx(k)*sys%m%h(1) * hzpsi(k, 1:sys%st%dim)
!!$               xzpsi(k, 1:sys%st%dim, 2) = sys%m%ly(k)*sys%m%h(2) * hzpsi(k, 1:sys%st%dim)
!!$               xzpsi(k, 1:sys%st%dim, 3) = sys%m%lz(k)*sys%m%h(3) * hzpsi(k, 1:sys%st%dim)
!!$            enddo
!!$
!!$            do j = 1, 3
!!$               vnl_xzpsi = (0.0_r8, 0.0_r8)
!!$               call zvnlpsi( sys, xzpsi(0:sys%m%np, 1:sys%st%dim, j), vnl_xzpsi(1:sys%m%np, 1:sys%st%dim))
!!$               do idim = 1, sys%st%dim
!!$                  x(j) = x(j) + 2*sys%st%occ(ist, ik)* zmesh_dotp(sys%m, R_CONJ(sys%st%zpsi(1:sys%m%np, idim, ist, ik)), &
!!$                                            vnl_xzpsi(1:sys%m%np, idim) )
!!$               enddo
!!$            enddo
!!$
!!$        enddo
!!$    enddo
!!$    deallocate(xzpsi, hzpsi, hhzpsi, vnl_xzpsi)
!!$
!!$#if defined(HAVE_MPI) && defined(MPI_TD)
!!$    if(present(reduce) .and. reduce) then
!!$       call MPI_ALLREDUCE(x(1), y(1), 3, &
!!$                          MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
!!$       x = y
!!$    end if
!!$#endif
!!$    acc = x
!!$
!!$    call pop_sub(); return
!!$  end subroutine td_calc_tacc2
