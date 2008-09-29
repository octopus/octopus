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
! This routine fills state psi with an atomic orbital -- provided
! by the pseudopotential structure in geo.
! ---------------------------------------------------------
subroutine X(lcao_atomic_orbital) (this, iorb, m, h, geo, sb, psi, spin_channel)
  type(lcao_t),             intent(in)    :: this
  integer,                  intent(in)    :: iorb
  type(mesh_t),             intent(in)    :: m
  type(simul_box_t),        intent(in)    :: sb
  type(hamiltonian_t),      intent(in)    :: h
  type(geometry_t), target, intent(in)    :: geo
  R_TYPE,                   intent(inout) :: psi(:, :)
  integer,                  intent(in)    :: spin_channel

  type(species_t), pointer :: s
  type(periodic_copy_t)   :: pc
  integer :: icell, idim, iatom, jj, ip
  FLOAT :: pos(MAX_DIM)
  FLOAT, allocatable :: ao(:)
  type(profile_t), save :: prof

  call profiling_in(prof, "ATOMIC_ORBITAL")
  call push_sub('lcao_inc.Xlcao_atomic_orbital')

  ASSERT(iorb >= 1)
  ASSERT(iorb <= this%maxorbs)

  psi(1:m%np, 1:h%d%dim) = R_TOTYPE(M_ZERO)

  iatom = this%atom(iorb)
  jj = this%level(iorb)
  idim = this%ddim(iorb)
  s => geo%atom(iatom)%spec
  ASSERT(jj <= s%niwfs)

  ALLOCATE(ao(1:m%np), m%np)

  if (.not. simul_box_is_periodic(sb)) then

    call species_get_orbital(s, m, jj, calc_dim, max(spin_channel, idim), geo%atom(iatom)%x, ao)

    do ip = 1, m%np
      psi(ip, idim) = ao(ip)
    end do

  else

    call periodic_copy_init(pc, sb, geo%atom(iatom)%x, range = species_get_iwf_radius(s, jj, spin_channel))
    do icell = 1, periodic_copy_num(pc)
      pos = periodic_copy_position(pc, sb, icell)

      call species_get_orbital(s, m, jj, calc_dim, max(spin_channel, idim), pos, ao)
      
      do ip = 1, m%np
        psi(ip, idim) = psi(ip, idim) + ao(ip)
      end do

    end do
    call periodic_copy_end(pc)

  end if

  call pop_sub()
  call profiling_out(prof)

end subroutine X(lcao_atomic_orbital)

! ---------------------------------------------------------
subroutine X(lcao_wf) (this, st, gr, geo, h, start)
  type(lcao_t),        intent(inout) :: this
  type(states_t),      intent(inout) :: st
  type(grid_t),        intent(inout) :: gr
  type(geometry_t),    intent(in)    :: geo
  type(hamiltonian_t), intent(in)    :: h
  integer,             intent(in)    :: start

  integer :: nst, ik, n1, n2, idim, lcao_start
  R_TYPE, allocatable :: hpsi(:, :, :), overlap(:, :, :)
  FLOAT, allocatable :: ev(:)
  R_TYPE, allocatable :: hamilt(:, :, :), lcaopsi(:, :, :), lcaopsi2(:, :)
  integer :: kstart, kend, ispin
  
  call push_sub('lcao_inc.Xlcao_wf')
  
  nst = st%nst
  kstart = st%d%kpt%start
  kend = st%d%kpt%end

  ! Allocation of variables

  ALLOCATE(lcaopsi(1:NP, 1:st%d%dim, 1:st%d%spin_channels), NP*st%d%dim*st%d%spin_channels)
  ALLOCATE(lcaopsi2(1:NP, 1:st%d%dim), NP*st%d%dim)
  ALLOCATE(hpsi(NP, st%d%dim, kstart:kend), NP*st%d%dim*st%d%kpt%nlocal)
  ALLOCATE(hamilt(this%norbs, this%norbs, kstart:kend), this%norbs**2*st%d%kpt%nlocal)
  ALLOCATE(overlap(this%norbs, this%norbs, st%d%spin_channels), this%norbs**2*st%d%spin_channels)

  do n1 = 1, this%norbs
    
    do ispin = 1, st%d%spin_channels
      call X(lcao_atomic_orbital)(this, n1, gr%m, h, geo, gr%sb, lcaopsi(:, :, ispin), ispin)
    end do

    do ik = kstart, kend
      ispin = states_dim_get_spin_index(st%d, ik)
      call X(hpsi)(h, gr, lcaopsi(:, :, ispin), hpsi(:, :, ik), n1, ik)
    end do

    do n2 = n1, this%norbs
      do ispin = 1, st%d%spin_channels
        call X(lcao_atomic_orbital)(this, n2, gr%m, h, geo, gr%sb, lcaopsi2, ispin)
        overlap(n1, n2, ispin) = X(mf_dotp)(gr%m, st%d%dim, lcaopsi(:, :, ispin), lcaopsi2)
        overlap(n2, n1, ispin) = R_CONJ(overlap(n1, n2, ispin))
        do ik = kstart, kend
          if(ispin /= states_dim_get_spin_index(st%d, ik)) cycle
          hamilt(n1, n2, ik) = X(mf_dotp)(gr%m, st%d%dim, hpsi(:, :, ik), lcaopsi2)
          hamilt(n2, n1, ik) = R_CONJ(hamilt(n1, n2, ik))
        end do
      end do
    end do
  end do

  deallocate(hpsi)

  ALLOCATE(ev(this%norbs), this%norbs)

  lcao_start = start
  if(st%parallel_in_states .and. st%st_start > start) lcao_start = st%st_start

  do ik =  kstart, kend
    ispin = states_dim_get_spin_index(st%d, ik)
    call lalg_geneigensolve(this%norbs, hamilt(1:this%norbs, 1:this%norbs, ik), overlap(:, :, ispin), ev)

    st%eigenval(start:nst, ik) = ev(start:nst)

    st%X(psi)(1:NP, 1:st%d%dim, lcao_start:st%st_end, ik) = R_TOTYPE(M_ZERO)
 
    ! Change of base
    do n2 = 1, this%norbs
      call X(lcao_atomic_orbital)(this, n2, gr%m, h, geo, gr%sb, lcaopsi2, ispin)
      
      do idim = 1, st%d%dim
        do n1 = lcao_start, st%st_end
          call lalg_axpy(NP, hamilt(n2, n1, ik), lcaopsi2(:, idim), st%X(psi)(:, idim, n1, ik))
        end do
      end do
      
    end do
    
  end do

  deallocate(ev, hamilt)
  
  call pop_sub()
end subroutine X(lcao_wf)

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
