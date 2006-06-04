!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch, X. Andrade.
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
subroutine X(lr_output) (st, gr, lr, dir, tag, outp)
  type(states_t),   intent(inout) :: st
  type(grid_t),     intent(inout) :: gr
  type(lr_t),       intent(inout) :: lr
  character(len=*), intent(in)    :: dir
  integer,          intent(in)    :: tag
  type(output_t),   intent(in)    :: outp

  integer :: ik, ist, idim, ierr
  character(len=80) :: fname
  FLOAT :: u
  FLOAT, allocatable :: dtmp(:)
  

  call push_sub('lr_response_out.lr_output')

  u = M_ONE/units_out%length%factor**NDIM

  if(iand(outp%what, output_density).ne.0) then
    do ist = 1, st%d%nspin
      write(fname, '(a,i1,a,i1)') 'lr_density-', ist, '-', tag
      call X(output_function)(outp%how, dir, fname, gr%m, gr%sb, lr%X(dl_rho)(:, ist), u, ierr)
    end do
  end if

!#ifdef COMPLEX_WFNS
!  if(iand(outp%what, output_current).ne.0) then
!    ! calculate current first
!    call calc_paramagnetic_current(gr, st, st%j)
!    do is = 1, st%d%nspin
!      do idim = 1, NDIM
!        write(fname, '(a,i1,a,a)') 'current-', is, '-', index2axis(idim)
!        call doutput_function(outp%how, dir, fname, gr%m, gr%sb, st%j(:, idim, is), u, ierr)
!      end do
!    end do
!  end if
!#endif

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
    
    integer :: i, is, ist, idim
    R_TYPE, allocatable :: gpsi(:,:), gdl_psi(:,:)
    FLOAT,  allocatable :: rho(:), grho(:,:), dl_rho(:), gdl_rho(:,:)
    FLOAT,  allocatable :: dl_elf(:,:), elf(:,:), de(:,:)
    FLOAT :: dl_d0, d0
    FLOAT :: f, s

    FLOAT, parameter :: dmin = CNST(1e-10)
    FLOAT :: u, ik_weight

    ALLOCATE(   gpsi(NP, NDIM), NP*NDIM)
    ALLOCATE(gdl_psi(NP, NDIM), NP*NDIM)
    ALLOCATE(   grho(NP, NDIM), NP*NDIM)
    ALLOCATE(gdl_rho(NP, NDIM), NP*NDIM)

    ALLOCATE(   rho(NP), NP)
    ALLOCATE(dl_rho(NP), NP)

    ALLOCATE(   elf(NP, st%d%nspin), NP*st%d%nspin)
    ALLOCATE(    de(NP, st%d%nspin), NP*st%d%nspin)
    ALLOCATE(dl_elf(NP, st%d%nspin), NP*st%d%nspin)
    
    !calculate the gs elf
    call states_calc_elf(st, gr, elf, de)

    ! single or double occupancy
    if(st%d%nspin == 1) then
      s = M_TWO
    else
      s = M_ONE
    end if

    dl_elf = M_ZERO

    do is = 1, st%d%nspin
      rho  = M_ZERO
      grho = M_ZERO
      
      dl_rho  = M_ZERO
      gdl_rho = M_ZERO

      do ik = is, st%d%nik, st%d%nspin
        do ist = 1, st%nst
          ik_weight = st%d%kweights(ik)*st%occ(ist, ik)/s

          do idim = 1, st%d%dim
          
            call X(f_gradient)(gr%sb, gr%f_der, st%X(psi)   (:, idim, ist, is), gpsi)
            call X(f_gradient)(gr%sb, gr%f_der, lr%X(dl_psi)(:, idim, ist, is), gdl_psi)

            ! sum over states to obtain the spin-density
            rho(1:NP)    = rho(1:NP)    + ik_weight * abs(st%X(psi)(1:NP, idim, ist, is))**2
            dl_rho(1:NP) = dl_rho(1:NP) + ik_weight * M_TWO *  &
               R_REAL(R_CONJ(st%X(psi)(1:NP, idim, ist, is)) * lr%X(dl_psi)(1:NP, idim, ist, is))

            do i = 1, NDIM
              grho(1:NP,i)    = grho(1:NP, i)   + ik_weight *  &
                 M_TWO * R_REAL(R_CONJ(st%X(psi)(1:NP, idim, ist, is))*gpsi(1:NP,i))
              
              gdl_rho(1:NP,i) = gdl_rho(1:NP,i) + ik_weight * M_TWO * ( &
                 R_REAL(R_CONJ(st%X(psi)(1:NP, idim, ist, is)) * gdl_psi(1:NP,i)) +      &
                 R_REAL(gpsi(1:NP,i) * R_CONJ(lr%X(dl_psi)(1:NP, idim, ist, is)))  )
                 
                 
              !jj(:,i)   =  jj(:,i)  + ik_weight *  &
              !   aimag(conjg(wf_psi(:))*gwf_psi(:,i))
            end do

            do i = 1, NP
              dl_elf(i, is) = dl_elf(i, is) +                             &
                 dl_rho(i) * ik_weight * sum(R_ABS(gpsi(i, 1:NDIM))**2) + &
                 rho(i)    * ik_weight * M_TWO*R_REAL(sum(R_CONJ(gpsi(i,1:NDIM))*gdl_psi(i,1:NDIM)))
            end do
          
          end do
        end do
      end do
      
      do i= 1, NP
        if(abs(st%rho(i, is)) >= dmin) then
          dl_elf(i, is) = dl_elf(i, is)                       &
             - M_HALF * sum(grho(i, 1:NDIM)*gdl_rho(i, 1:NDIM))
        end if
      end do
      
      write(fname, '(a,a,i1,a,i1)') trim(filename1), '-', is, '-', tag 
      call doutput_function(outp%how, dir, trim(fname), gr%m, gr%sb, dl_elf(:,is), M_ONE, ierr)
      
      !normalization 
      f = M_THREE/M_FIVE*(M_SIX*M_PI**2)**M_TWOTHIRD
      do i = 1, NP

        if(abs(st%rho(i, is)) >= dmin) then
          d0    = f * rho(i)**(M_EIGHT/M_THREE) + CNST(1e-5)
          dl_d0 = M_EIGHT/M_THREE * f * dl_rho(i) * rho(i)**(M_FIVE/M_THREE)
        
          dl_elf(i, is) = - M_TWO * elf(i,is)**2 * (de(i, is)/d0) * &
             ((dl_elf(i, is) - (de(i, is)/d0)*dl_d0)/d0)
        else
          dl_elf(i, is) = M_ZERO
        end if

      end do
        
      write(fname, '(a,a,i1,a,i1)') trim(filename2), '-', is, '-', tag 
      call doutput_function(outp%how, dir, trim(fname), gr%m, gr%sb, dl_elf(:,is), M_ONE, ierr)
      
    end do

    deallocate(gpsi)
    deallocate(gdl_psi)
    deallocate(grho)

    deallocate(elf)
    deallocate(de)
    deallocate(dl_elf)
    
  end subroutine lr_elf

end subroutine X(lr_output)

