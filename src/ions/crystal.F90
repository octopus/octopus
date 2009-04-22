!! $Id$


#include "global.h"

module crystal_m
  use blas_m
  use global_m
  use geometry_m
  use lalg_adv_m
  use messages_m
  use simul_box_m
  use symmetries_m
  use profiling_m
  use units_m

  implicit none

  private
  
  public :: crystal_init
  
contains

  subroutine crystal_monkhorstpack_generate(sb, naxis, shift, nkpoints, kpoints)
    type(simul_box_t), intent(in)  :: sb
    integer,           intent(in)  :: naxis(1:MAX_DIM)
    FLOAT,             intent(in)  :: shift(1:MAX_DIM)
    integer,           intent(out) :: nkpoints
    FLOAT,             intent(out) :: kpoints(:, :)

    ! Implements the Monkhorst-Pack scheme.

    ! Sets up uniform array of k points. Use the normal MP scheme
    ! (PRB13, 5188, (1976)) when sx=sy=sz=0.5. If sx=sy=0,
    ! the special hexagonal scheme is used (PRB16, 1748, (1977))

    FLOAT :: dx(1:MAX_DIM)
    integer :: ii, jj, kk, idir, ix(1:MAX_DIM)

    call push_sub('crystal.crystal_monkhorstpack_generate')

    ! nx, ny, and nz are the number of points in the three
    ! directions dermined by the lattice wave vectors. sx, sy, and
    ! sz shift the grid of integration points from the origin.

    !generate k-points using the MP scheme
    
    dx(1:MAX_DIM) = M_ONE/real(2*naxis(1:MAX_DIM), REAL_PRECISION)

    kpoints = M_ZERO
    nkpoints = 0
    ix=1
    do ii = 1, naxis(1)
      do jj = 1, naxis(2)
        do kk = 1, naxis(3)
          ix(1:3) = (/ii, jj, kk/)
          nkpoints = nkpoints + 1
          forall(idir = 1:sb%periodic_dim)
            kpoints(idir, nkpoints) = (real(2*ix(idir) - naxis(idir) - 1, REAL_PRECISION) + shift(idir))*dx(idir)
          end forall
        end do
      end do
    end do
    
    call pop_sub()

  end subroutine crystal_monkhorstpack_generate

  subroutine crystal_init(geo, sb, nk_axis, shift, use_symmetries, use_time_reversal, nkpoints, kpoints, weights)
    type(geometry_t),  intent(in)    :: geo
    type(simul_box_t), intent(in)    :: sb
    integer,           intent(in)    :: nk_axis(MAX_DIM)
    FLOAT,             intent(in)    :: shift(MAX_DIM)
    logical,           intent(in)    :: use_symmetries
    logical,           intent(in)    :: use_time_reversal
    integer,           intent(inout) :: nkpoints
    FLOAT,             intent(out)   :: kpoints(:, :)
    FLOAT,             intent(out)   :: weights(:)

    integer :: ik
    type(symmetries_t) :: symm

    call push_sub('crystal.crystal_init')

    call crystal_monkhorstpack_generate(sb, nk_axis, shift, nkpoints, kpoints)
    
    if(use_symmetries) then
      call symmetries_init(symm, geo, sb)
      call crystal_monkhorstpack_reduce(symm, use_time_reversal, nkpoints, kpoints, weights)
      call symmetries_end(symm)
    else
      forall(ik = 1:nkpoints) weights(ik) = M_ONE/dble(nkpoints)
    end if

    write(message(1),'(a)')
    write(message(2),'(1x,i3,a)') nkpoints, ' k points generated from parameters :'
    write(message(3),'(1x,a)') '---------------------------------------------------'
    write(message(4),'(4x,a,3i5,6x,a,3f6.2)') 'n =', nk_axis(1:3), 's = ', shift(1:3)
    write(message(5),'(a)')
    call write_info(5)

    do ik = 1, nkpoints
      write(message(1), '(i5,4f12.6)') ik, weights(ik), kpoints(1:3, ik)
    call write_info(1)
    end do

    call pop_sub()
  end subroutine crystal_init
  
  subroutine crystal_monkhorstpack_reduce(symm, time_reversal, nkpoints, kpoints, weights)
    type(symmetries_t), intent(in)    :: symm
    logical,            intent(in)    :: time_reversal
    integer,            intent(inout) :: nkpoints
    FLOAT,              intent(inout) :: kpoints(:, :)
    FLOAT,              intent(out)   :: weights(:)

    integer :: nreduced
    FLOAT, allocatable :: reduced(:, :)
    
    FLOAT :: dw
    integer ik, iop, ik2
    FLOAT :: tran(MAX_DIM), tran_inv(MAX_DIM)
    integer, allocatable :: kmap(:)

    call push_sub('crystal.crystal_monkhortstpack_reduce')

    ! reduce to irreducible zone

    ! kmap is used to mark reducible k points and also to
    ! map reducible to irreducible k points

    SAFE_ALLOCATE(kmap(1:nkpoints))
    SAFE_ALLOCATE(reduced(1:MAX_DIM, 1:nkpoints))

    forall(ik = 1:nkpoints) kmap(ik) = ik
    
    dw = M_ONE/real(nkpoints, REAL_PRECISION)

    nreduced = 0

    do ik = 1, nkpoints
      if (kmap(ik) /= ik) cycle
      
      ! new irreducible point
      ! mareduced with negative kmap
      
      nreduced = nreduced + 1
      reduced(1:MAX_DIM, nreduced) = kpoints(1:MAX_DIM, ik)
      
      kmap(ik) = -nreduced
      
      weights(nreduced) = dw

      if (ik == nkpoints) cycle
      
      ! operate with the symmetry operations
      
      do iop = 1, symmetries_number(symm)
        call symmetries_apply(symm, iop, reduced(:, nreduced), tran)
        tran_inv(1:MAX_DIM) = -tran(1:MAX_DIM)
           
        ! remove (mark) k points related to irreducible reduced by symmetry
        do ik2 = ik + 1, nkpoints
          if (kmap(ik2) /= ik2) cycle
          
          ! both the transformed rk ...
          if (all( abs(tran(1:3) - kpoints(1:3, ik2)) <= CNST(1.0e-5))) then
            weights(nreduced) = weights(nreduced) + dw
            kmap(ik2) = nreduced
            cycle
          end if

          ! and its inverse
          if (time_reversal .and. all(abs(tran_inv(1:3) - kpoints(1:3, ik2)) <= CNST(1.0e-5)) ) then
            weights(nreduced) = weights(nreduced) + dw
            kmap(ik2) = nreduced
          end if
          
        end do
      end do
    end do
    
    nkpoints = nreduced
    do ik = 1, nreduced
      kpoints(1:3, ik) = reduced(1:3, ik)
    end do

    call pop_sub()
  end subroutine crystal_monkhorstpack_reduce

end module crystal_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
