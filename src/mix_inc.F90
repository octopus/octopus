! Calculates the new density out the wavefunctions and occupations...
subroutine R_FUNC(calcdens)(st, np, rho, reduce)
  type(states_type), intent(IN) :: st
  integer, intent(in) :: np
  real(r8), intent(out) :: rho(np, st%ispin)
  logical, intent(in), optional :: reduce

  integer :: i, ik, p, sp
#ifdef HAVE_MPI
  real(r8), allocatable :: reduce_rho(:,:)
  integer :: ierr
#endif

  if(st%ispin == 2) then
    sp = 2
  else
    sp = 1
  end if

  ! TODO: for polymers, do not forget to introduce integration factor
  ! in momentum space
  rho = 0._r8
  do ik = 1, st%nik, sp
    do p  = st%st_start, st%st_end
      do i = 1, np
        ! spin-up density
        rho(i, 1) = rho(i, 1) + st%occ(p, ik)*R_ABS(st%R_FUNC(psi)(i, 1, p, ik))**2

        ! spin-down density
        if(st%ispin == 2) then
          rho(i, 2) = rho(i, 2) + st%occ(p, ik+1)*R_ABS(st%R_FUNC(psi)(i, 1, p, ik+1))**2
        end if

        ! off-diagonal densities
        if(st%ispin == 4) then
          rho(i, 2) = rho(i, 2) + st%occ(p, ik)*R_ABS(st%R_FUNC(psi)(i, 2, p, ik))**2
!          rho(i, 3) = st%occ(p, ik)*&
!               st%R_FUNC(psi)(i, 1, p, ik)*R_CONJ(st%R_FUNC(psi)(i, 2, p, ik))
!          rho(i, 4) = R_CONJ(rho(i, 3)) ! this is in principle not necessary!
        end if

      end do
    end do
  end do

#if defined(HAVE_MPI) && defined(MPI_TD)
  ! reduce density (assumes memory is contiguous)
  if(present(reduce) .and. reduce) then
    i = min(2, st%ispin)
    allocate(reduce_rho(1:m%np, i))
    call MPI_ALLREDUCE(rho(1, 1), reduce_rho(1, 1), m%np*i, &
         MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    rho = reduce_rho
    ! TODO update also off diagonal...
    deallocate(reduce_rho)
  end if
#endif

  return
end subroutine R_FUNC(calcdens)
