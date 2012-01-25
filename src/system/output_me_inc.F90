!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!!
!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any %later version.
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
!! $Id: unocc.F90 5320 2009-04-23 15:29:57Z marques $


! ---------------------------------------------------------
!> Prints out the multipole matrix elements between KS states.
!!
!! It prints the states to the file opened in iunit.
!! It prints the (ll,mm) multipole moment, for
!! the Kohn-Sham states in the irreducible subspace ik.
! ---------------------------------------------------------
subroutine X(output_me_ks_multipoles)(fname, st, gr, ll, mm, ik)
  character(len=*), intent(in) :: fname
  type(states_t),   intent(in) :: st
  type(grid_t),     intent(in) :: gr
  integer,          intent(in) :: ll, mm, ik
  
  integer :: ist, jst, ip, iunit
  FLOAT, allocatable :: multipole(:)
  FLOAT :: rr, xx(MAX_DIM), ylm
  R_TYPE :: multip_element
    R_TYPE, allocatable :: psii(:, :), psij(:, :)

  PUSH_SUB(X(output_me_ks_multipoles))

  SAFE_ALLOCATE(psii(1:gr%mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(psij(1:gr%mesh%np, 1:st%d%dim))
  
  iunit = io_open(file = fname, action = 'write')

  write(iunit, fmt = '(a)') '# Multipole matrix elements file: <Phi_i | r**l * Y_{lm}(theta,phi) | Phi_j>' 
  write(iunit, fmt = '(a,i2,a,i2)') '# l =', ll, '; m =', mm
  write(iunit, fmt = '(a,i4)')      '# ik =', ik
  if(ll>1) then
    write(iunit, fmt = '(a,i1)') '# Units = ['//trim(units_abbrev(units_out%length))//']^', ll
  else
    write(iunit, fmt = '(a)')    '# Units = ['//trim(units_abbrev(units_out%length))//']'
  end if
  
  SAFE_ALLOCATE(multipole(1:gr%mesh%np_part))

  multipole = M_ZERO
  do ip = 1, gr%mesh%np
    call mesh_r(gr%mesh, ip, rr, coords = xx)
    call loct_ylm(1, xx(1), xx(2), xx(3), ll, mm, ylm)
    multipole(ip) = rr**ll * ylm
  end do
  
  do ist = 1, st%nst

    call states_get_state(st, gr%mesh, ist, 1, psii)

    do jst = 1, st%nst

      call states_get_state(st, gr%mesh, jst, 1, psij)

      psij(1:gr%mesh%np, 1) = psij(1:gr%mesh%np, 1)*multipole(1:gr%mesh%np)

      multip_element = X(mf_dotp)(gr%mesh, st%d%dim, psii, psij)

      multip_element = units_from_atomic(units_out%length**ll, multip_element)

      write(iunit, fmt='(f20.12)', advance = 'no') R_REAL(multip_element)
#if defined(R_TCOMPLEX)
      write(iunit, fmt='(a1,f20.12)', advance = 'no') ',', R_AIMAG(multip_element)
#endif
      write(iunit, fmt='(a)', advance = 'no') '   '
    end do
    write(iunit, '(a)') ''
  end do

  SAFE_DEALLOCATE_A(psii)
  SAFE_DEALLOCATE_A(psij)
  SAFE_DEALLOCATE_A(multipole)
  call io_close(iunit)

  POP_SUB(X(output_me_ks_multipoles))
end subroutine X(output_me_ks_multipoles)


! ---------------------------------------------------------
subroutine X(one_body) (dir, gr, geo, st, hm)
  character(len=*),    intent(in)    :: dir
  type(grid_t),        intent(inout) :: gr
  type(geometry_t),    intent(in)    :: geo
  type(states_t),      intent(inout) :: st
  type(hamiltonian_t), intent(in)    :: hm

  integer ist, jst, iunit, idir, iatom, np
  R_TYPE :: me, exp_r, exp_g, corr
  R_TYPE, allocatable :: gpsi(:,:), cpsi(:,:)
  R_TYPE, allocatable :: psii(:, :), psij(:, :)

  SAFE_ALLOCATE(psii(1:gr%mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(psij(1:gr%mesh%np, 1:st%d%dim))
  
  PUSH_SUB(X(one_body))

  np = gr%mesh%np

  iunit = io_open(trim(dir)//'/output_me_one_body', action='write')

  do ist = 1, st%nst

    call states_get_state(st, gr%mesh, ist, 1, psii)
    
    do jst = 1, st%nst
      if(jst > ist) cycle
      
      call states_get_state(st, gr%mesh, jst, 1, psij)

      psij(1:np, 1) = R_CONJ(psii(1:np, 1))*hm%Vhxc(1:np, 1)*psij(1:np, 1)

      me = st%eigenval(ist,1) - X(mf_integrate)(gr%mesh, psij(:, 1))

      write(iunit, *) ist, jst, me
    end do
  end do

  SAFE_ALLOCATE(gpsi(1:gr%mesh%np_part, 1:gr%sb%dim))
  SAFE_ALLOCATE(cpsi(1:gr%mesh%np_part, 1:1))

  iunit = io_open(trim(dir)//'/output_me_gauge', action='write')

  do ist = 1, st%nst

    call states_get_state(st, gr%mesh, ist, 1, psii)

    do jst = 1, st%nst
      if(st%occ(ist, 1) < CNST(0.0001)) cycle
      if(st%occ(jst, 1) > CNST(0.0001)) cycle

      call states_get_state(st, gr%mesh, jst, 1, psij)

      call X(derivatives_grad)(gr%der, psij(:, 1), gpsi)
       
      do idir = 1, 3
         exp_r = X(mf_integrate)(gr%mesh, R_CONJ(psii(1:np, 1))*gr%mesh%x(1:np, idir)*psij(1:np, 1))
         exp_g = X(mf_integrate)(gr%mesh, R_CONJ(psii(1:np, 1))*gpsi(1:np, idir))
         
         corr = M_ZERO
         do iatom = 1, geo%natoms
           cpsi = M_ZERO
           call X(projector_commute_r)(hm%ep%proj(iatom), gr, 1, idir, 1, psij, cpsi)
           corr = corr + X(mf_integrate)(gr%mesh, R_CONJ(psii(1:np, 1))*cpsi(1:np, 1))
         end do

         me = (st%eigenval(jst, 1) - st%eigenval(ist, 1))*exp_r
         
         write(iunit, *) ist, jst, idir, me, me - (exp_g + corr)

       end do
       
     end do
   end do
   
  SAFE_DEALLOCATE_A(gpsi)
  SAFE_DEALLOCATE_A(psii)
  SAFE_DEALLOCATE_A(psij)

  call io_close(iunit)
  POP_SUB(X(one_body))
end subroutine X(one_body)


! ---------------------------------------------------------
subroutine X(two_body) (dir, gr, st)
  character(len=*), intent(in)    :: dir
  type(grid_t),     intent(inout) :: gr
  type(states_t),   intent(in)    :: st

  integer ist, jst, kst, lst, iunit
  R_TYPE :: me
  R_TYPE, allocatable :: nn(:), vv(:)
  R_TYPE, allocatable :: psii(:, :), psij(:, :), psik(:, :), psil(:, :)

  PUSH_SUB(X(two_body))

  iunit = io_open(trim(dir)//'/output_me_two_body', action='write')

  SAFE_ALLOCATE(nn(1:gr%mesh%np))
  SAFE_ALLOCATE(vv(1:gr%mesh%np))
  SAFE_ALLOCATE(psii(1:gr%mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(psij(1:gr%mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(psik(1:gr%mesh%np, 1:st%d%dim))
  SAFE_ALLOCATE(psil(1:gr%mesh%np, 1:st%d%dim))

  do ist = 1, st%nst

    call states_get_state(st, gr%mesh, ist, 1, psii)

    do jst = 1, st%nst
      if(jst > ist) cycle

      call states_get_state(st, gr%mesh, jst, 1, psij)

      nn(1:gr%mesh%np) = R_CONJ(psii(1:gr%mesh%np, 1))*psij(1:gr%mesh%np, 1)
      call X(poisson_solve)(psolver, vv, nn, all_nodes=.false.)

      do kst = 1, st%nst
        if(kst > ist) cycle
        
        call states_get_state(st, gr%mesh, kst, 1, psik)

        do lst = 1, st%nst
          if(lst > kst) cycle
          if(lst > jst) cycle

          call states_get_state(st, gr%mesh, lst, 1, psil)

          psil(1:gr%mesh%np, 1) = vv(1:gr%mesh%np)*psik(1:gr%mesh%np, 1)*R_CONJ(psil(1:gr%mesh%np, 1))

          me = X(mf_integrate)(gr%mesh, psil(:, 1))

          write(iunit, *) ist, jst, kst, lst, me
        end do
      end do
    end do
  end do

  SAFE_DEALLOCATE_A(psii)
  SAFE_DEALLOCATE_A(psij)
  SAFE_DEALLOCATE_A(psik)
  SAFE_DEALLOCATE_A(psil)

  call io_close(iunit)
  POP_SUB(X(two_body))
end subroutine X(two_body)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
