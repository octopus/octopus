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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!


! ---------------------------------------------------------
!> calculates the eigenvalues of the orbitals
subroutine X(calculate_eigenvalues)(hm, der, st, time)
  type(hamiltonian_t), intent(in)    :: hm
  type(derivatives_t), intent(inout) :: der
  type(states_t),      intent(inout) :: st
  FLOAT,   optional,   intent(in)    :: time

  R_TYPE, allocatable :: eigen(:, :)
  logical :: cmplxscl

  PUSH_SUB(X(calculate_eigenvalues))
  
  cmplxscl = hm%cmplxscl%space

  if(hm%theory_level == CLASSICAL) then
    st%eigenval = M_ZERO
    POP_SUB(X(calculate_eigenvalues))
    return
  end if

  if(debug%info) then
    write(message(1), '(a)') 'Debug: Calculating eigenvalues.'
    call messages_info(1)
  end if

  st%eigenval = M_ZERO
  if(cmplxscl) st%zeigenval%Im = M_ZERO

  SAFE_ALLOCATE(eigen(st%st_start:st%st_end, st%d%kpt%start:st%d%kpt%end))
  call X(calculate_expectation_values)(hm, der, st, eigen, time = time)

  st%eigenval(st%st_start:st%st_end, st%d%kpt%start:st%d%kpt%end) = &
    real(eigen(st%st_start:st%st_end, st%d%kpt%start:st%d%kpt%end), REAL_PRECISION)
#ifdef R_TCOMPLEX    
  if(cmplxscl) st%zeigenval%Im(st%st_start:st%st_end, st%d%kpt%start:st%d%kpt%end) = &
    aimag(eigen(st%st_start:st%st_end, st%d%kpt%start:st%d%kpt%end))
#endif

  call comm_allreduce(st%st_kpt_mpi_grp%comm, st%eigenval)
  if(cmplxscl) call comm_allreduce(st%st_kpt_mpi_grp%comm, st%zeigenval%Im)

  SAFE_DEALLOCATE_A(eigen)

  POP_SUB(X(calculate_eigenvalues))
end subroutine X(calculate_eigenvalues)

subroutine X(calculate_expectation_values)(hm, der, st, eigen, time, terms)
  type(hamiltonian_t), intent(in)    :: hm
  type(derivatives_t), intent(inout) :: der
  type(states_t),      intent(inout) :: st
  R_TYPE,              intent(out)   :: eigen(st%st_start:, st%d%kpt%start:) !< (:st%st_end, :st%d%kpt%end)
  FLOAT,   optional,   intent(in)    :: time
  integer, optional,   intent(in)    :: terms

  integer :: ik, minst, maxst, ib
  type(batch_t) :: hpsib
  type(profile_t), save :: prof
  logical :: cmplxscl

  PUSH_SUB(X(calculate_expectation_values))
  
  call profiling_in(prof, "EIGENVALUE_CALC")

  cmplxscl = hm%cmplxscl%space

  do ik = st%d%kpt%start, st%d%kpt%end
    do ib = st%group%block_start, st%group%block_end

      minst = states_block_min(st, ib)
      maxst = states_block_max(st, ib)

      call batch_copy(st%group%psib(ib, ik), hpsib)

      if(hamiltonian_apply_packed(hm, der%mesh)) then
        call batch_pack(st%group%psib(ib, ik))
        if(st%have_left_states) call batch_pack(st%psibL(ib, ik))
        call batch_pack(hpsib, copy = .false.)
      end if

      call X(hamiltonian_apply_batch)(hm, der, st%group%psib(ib, ik), hpsib, ik, time = time, terms = terms)
      if(st%have_left_states) then
        call X(mesh_batch_dotp_vector)(der%mesh, st%psibL(ib, ik), hpsib, eigen(minst:maxst, ik), cproduct = cmplxscl)
      else
        call X(mesh_batch_dotp_vector)(der%mesh, st%group%psib(ib, ik), hpsib, eigen(minst:maxst, ik), cproduct = cmplxscl)        
      end if
      if(hamiltonian_apply_packed(hm, der%mesh)) then
        call batch_unpack(st%group%psib(ib, ik), copy = .false.)
        if(st%have_left_states) call batch_unpack(st%psibL(ib, ik), copy = .false.)
      end if

      call batch_end(hpsib, copy = .false.)

    end do
  end do

  call profiling_out(prof)
  POP_SUB(X(calculate_expectation_values))
end subroutine X(calculate_expectation_values)

! ---------------------------------------------------------
R_TYPE function X(energy_calc_electronic)(hm, der, st, terms) result(energy)
  type(hamiltonian_t), intent(in)    :: hm
  type(derivatives_t), intent(inout) :: der
  type(states_t),      intent(inout) :: st
  integer,             intent(in)    :: terms

  R_TYPE, allocatable  :: tt(:, :)
 
  PUSH_SUB(X(energy_calc_electronic))

  SAFE_ALLOCATE(tt(st%st_start:st%st_end, st%d%kpt%start:st%d%kpt%end))

  call X(calculate_expectation_values)(hm, der, st, tt, terms = terms)

  if(hm%cmplxscl%space) then
#ifdef R_TCOMPLEX
    energy = zstates_eigenvalues_sum(st, tt)
#else
    message(1) = "Internal error in energy_calc_electronic, real states but complex scaling."
    call messages_fatal(1)
#endif
  else
    energy = states_eigenvalues_sum(st, real(tt, REAL_PRECISION))
  end if

  SAFE_DEALLOCATE_A(tt)
  POP_SUB(X(energy_calc_electronic))
end function X(energy_calc_electronic)

!-----------------------------------------------

subroutine X(v_resp_calc)(st, gr, mu, d, vresp)
  integer, intent(in)            :: d					
  type(states_t), intent(in)     :: st                                                           
  FLOAT, intent(in)              :: mu	
  type(grid_t), intent(in)       :: gr
  FLOAT, pointer                 :: vresp(:,:)						
  FLOAT, pointer                 :: density(:,:), occ(:,:)          
  R_TYPE, pointer                :: psi(:,:)                        
  integer                        :: ip, ist, iq, ib                     
  FLOAT                          :: Kx         
  FLOAT                          :: weight

  !The Kx constant depends on dimentionality
  !write(*,*) "Dimentionality: ", d
  select case(d) 
    case(3); Kx = CNST(0.382106112)
    case(2); Kx = CNST(0.450158158)                    
    case(1); call messages_not_implemented('GLLB in 1D')
  end select                                                                
  
  !Pointing some useful quantities                                  
  ! in order to simplify the notation                               
  occ => st%occ
  density => st%rho

  SAFE_ALLOCATE(vresp(gr%mesh%np, 1))
  SAFE_ALLOCATE(psi(gr%mesh%np, 1))
  vresp = M_ZERO
  do iq = st%d%kpt%start, st%d%kpt%end                !Sum over k-points
    do ib = st%group%block_start, st%group%block_end  !Sum over block states
      do ist = 1, st%group%psib(ib, iq)%nst           !Sum over the states of the block
        if(st%eigenval(st%group%psib(ib,iq)%states(ist)%ist,1)<mu) then !We do not calculate weight if mu < epsilon(ist,1)
          psi => st%group%psib(ib, iq)%states(ist)%X(psi)         !Point to the ist state of the ib block
          !write(*,*) "mu =", mu
          !write(*,*) "eigenval =", st%eigenval
          weight = Kx*st%d%kweights(iq)*occ(st%group%psib(ib,iq)%states(ist)%ist, 1)&
                        *sqrt(mu-st%eigenval(st%group%psib(ib,iq)%states(ist)%ist,1))
          write(*,*) weight
          forall(ip = 1:gr%mesh%np)  
            vresp(ip,1) = vresp(ip,1) + weight*(R_REAL(psi(ip, 1))**2+& 
                          R_AIMAG(psi(ip, 1))**2)/(density(ip,1) + 1E-10)
          end forall
        end if
      end do
    end do
  end do

  !write(*,*) "vresp(100) =", vresp(100,1)
end subroutine X(v_resp_calc)




!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
