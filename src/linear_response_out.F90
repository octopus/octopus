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
      write(fname, '(a,i1,a,i1)') 'lr-density-', ist, '-',tag
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
            write(fname, '(a,i3.3,a,i3.3,a,i1,a,i1)') 'lr-wf-', ik, '-', ist, '-', idim,'-', tag
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
            write(fname, '(a,i3.3,a,i3.3,a,i1)') 'sqm-lr-wf-', ik, '-', ist, '-', idim
            dtmp = abs(lr%X(dl_psi) (:, idim, ist, ik))**2
            call doutput_function (outp%how, dir, fname, gr%m, gr%sb, dtmp, u, ierr)
          end do
        end do
      end if
    end do
    deallocate(dtmp)
  end if

  if(NDIM==3) then
    if(iand(outp%what, output_elf).ne.0)    call lr_elf('lr-elf')
  end if

  call pop_sub()
contains

  ! ---------------------------------------------------------
  subroutine lr_elf(filename)
    character(len=*), intent(in) :: filename
    
    integer :: spin, i, ist, idim
    R_TYPE, allocatable :: gpsi(:,:), gdl_psi(:,:)
    FLOAT,  allocatable :: grho(:,:), dl_rho(:),gdl_rho(:,:)
    FLOAT,  allocatable :: dl_elf(:,:)


    ALLOCATE(gpsi(1:NP, 1:NDIM), NP*NDIM)
    ALLOCATE(gdl_psi(1:NP, 1:NDIM), NP*NDIM)
    ALLOCATE(grho(1:NP, 1:NDIM), NP*NDIM)
    ALLOCATE(gdl_rho(1:NP, 1:NDIM), NP*NDIM)

    ALLOCATE(dl_rho(1:NP_PART), NP_PART*NDIM)

    ALLOCATE(dl_elf(1:NP, 1:st%d%nspin), NP)
    
    dl_elf = M_ZERO

    do spin = 1, st%d%nspin
      

      do ist = 1, st%nst
        do idim = 1, st%d%dim
          
          call X(f_gradient)(gr%sb, gr%f_der, st%X(psi)(:, idim, ist, spin), gpsi)
          call X(f_gradient)(gr%sb, gr%f_der, lr%X(dl_psi)(:, idim, ist, spin), gdl_psi)
          
          do i = 1, NP
            dl_elf(i,spin) = dl_elf(i,spin) + M_TWO*R_REAL(sum(R_CONJ(gpsi(i,1:NDIM))*gdl_psi(i,1:NDIM)))
          end do
            
        end do
      end do

      dl_rho(1:NP) = lr%X(dl_rho)(1:NP,spin)
      dl_rho(NP+1:NP_PART) = M_ZERO

      call df_gradient(gr%sb, gr%f_der, dl_rho(:), gdl_rho)
      call df_gradient(gr%sb, gr%f_der, st%rho(:,spin), grho)

      do i= 1, NP
        dl_elf(i,spin) = dl_elf(i,spin)                                             &
             - M_HALF   * sum(grho(i,1:NDIM)*R_REAL(gpsi(i,1:NDIM)))/st%rho(i,spin) & 
             + M_FOURTH * sum(grho(i,1:NDIM)**2)*dl_rho(i)/(st%rho(i,spin)**2)
      end do
      
      write(fname, '(a,a,i1,a,i1)') trim(filename), '-', spin, '-', tag 
      call doutput_function(outp%how, dir, trim(fname), gr%m, gr%sb, dl_elf(:,1), M_ONE, ierr)

    end do

    deallocate(gpsi)
    deallocate(gdl_psi)
    deallocate(grho)

    deallocate(dl_elf)
    
    
#if 0
    FLOAT :: f, d, s
    integer :: i, is, ik
    CMPLX, allocatable :: psi_fs(:), gpsi(:,:)
    FLOAT, allocatable :: c(:), r(:), gradr(:,:), j(:,:)
    type(X(cf_t)) :: cf_tmp

    FLOAT, parameter :: dmin = CNST(1e-10)

    ! single or double occupancy
    if(st%d%nspin == 1) then
      s = M_TWO
    else
      s = M_ONE
    end if

    ALLOCATE(c(NP), NP)

    do_is: do is = 1, st%d%nspin
      ALLOCATE(    r(NP),       NP)
      ALLOCATE(gradr(NP, NDIM), NP*NDIM)
      ALLOCATE(    j(NP, NDIM), NP*NDIM)
      r = M_ZERO; gradr = M_ZERO; j  = M_ZERO
      c = M_ZERO

      ALLOCATE(psi_fs(NP_PART),  NP_PART)
      ALLOCATE(gpsi  (NP, NDIM), NP*NDIM)
      do ik = is, st%d%nik, st%d%nspin
        do ist = 1, st%nst
          do idim = 1, st%d%dim

            psi_fs(:) = cmplx(st%X(psi)(:, idim, ist, ik), KIND=PRECISION)
            call zf_gradient(gr%sb, gr%f_der, psi_fs(:), gpsi)

            r(:) = r(:) + st%d%kweights(ik)*st%occ(ist, ik) * abs(psi_fs(:))**2
            do i = 1, NDIM
              gradr(:,i) = gradr(:,i) + st%d%kweights(ik)*st%occ(ist, ik) *  &
                M_TWO * real(conjg(psi_fs(:))*gpsi(:,i))
              j (:,i) =  j(:,i) + st%d%kweights(ik)*st%occ(ist, ik) *  &
                aimag(conjg(psi_fs(:))*gpsi(:,i))
            end do

            do i = 1, NP
              if(r(i) >= dmin) then
                c(i) = c(i) + st%d%kweights(ik)*st%occ(ist, ik)/s * &
                  sum(abs(gpsi(i, 1:NDIM))**2)
              end if
            end do
          end do
        end do
      end do
      deallocate(psi_fs, gpsi)

      do i = 1, NP
        if(r(i) >= dmin) then
          c(i) = c(i) - (M_FOURTH*sum(gradr(i, 1:NDIM)**2) + sum(j(i, 1:NDIM)**2))/(s*r(i))
        end if
      end do

      ! normalization
      f = M_THREE/M_FIVE*(M_SIX*M_PI**2)**M_TWOTHIRD
      do i = 1, NP
        if(abs(r(i)) >= dmin) then
          d    = f*(r(i)/s)**(M_FIVE/M_THREE)
          c(i) = M_ONE/(M_ONE + (c(i)/d)**2)
        else
          c(i) = M_ZERO
        end if
      end do

      deallocate(r, gradr, j)
      write(fname, '(a,a,i1)') trim(filename), '-', is
      call doutput_function(outp%how, dir, trim(fname), gr%m, gr%sb, c, M_ONE, ierr)

    end do do_is

    deallocate(c)
#endif 


  end subroutine lr_elf

end subroutine X(lr_output)

