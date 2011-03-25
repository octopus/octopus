


#include "global.h"

module permutations_m

  use global_m
  use loct_m
  use math_m
  use messages_m
  use profiling_m

  implicit none

  private

  public :: permutations_init, &
            permutations_write, &
            permutations_copy, &
            permutations_nullify, &
            permutations_end, &
            permutations_t

  type permutations_t 
   integer :: nn, npermutations, npairs
   integer, pointer :: allpermutations(:,:)
   integer, pointer :: permsign(:)

  end type permutations_t 



  contains

  subroutine permutations_init (nn, this)

    integer, intent(in)                 :: nn
    type(permutations_t), intent(inout) :: this

    integer :: i1, order, oldperm, iperm, newpos

    PUSH_SUB(permutations_init)

    this%nn = nn
    this%npermutations = factorial(nn)
    SAFE_ALLOCATE (this%allpermutations(1:max(1,nn),1:this%npermutations))
    SAFE_ALLOCATE (this%permsign(1:this%npermutations))

    this%allpermutations(:, :) = -999
    do i1 = 1, nn
      this%allpermutations(i1, 1) = i1
    end do
    this%permsign(1) = 1
   
    iperm = 1
    do order = 2, nn
      do oldperm = 1, factorial(order-1)
        do newpos = order-1, 1, -1
          iperm = iperm + 1
          this%allpermutations(1:newpos-1,     iperm) = this%allpermutations(1:newpos-1, oldperm)
          this%allpermutations(newpos,         iperm) = order
          this%allpermutations(newpos+1:order, iperm) = this%allpermutations(newpos:order-1, oldperm)
          this%allpermutations(order+1:nn,     iperm) = this%allpermutations(order+1:nn, oldperm)
 
          this%permsign(iperm) = this%permsign(oldperm) * (-1)**(order-newpos)
        end do
      end do
    end do

    POP_SUB(permutations_init)
  end subroutine permutations_init

  subroutine permutations_write (this)
    type(permutations_t), intent(inout) :: this
    integer :: iperm
    PUSH_SUB(permutations_write)
    
    do iperm = 1, this%npermutations
      write (message(1), '(a,I7,a,I7,a,10I7)') 'permutation ', iperm, &
              ' sign ', this%permsign(iperm), '= ', this%allpermutations(:,iperm)
      call messages_info(1)
    end do
    
    POP_SUB(permutations_write)
  end subroutine permutations_write

  subroutine permutations_nullify (this)
    type(permutations_t), intent(inout) :: this
    PUSH_SUB(permutations_nullify)
    nullify(this%allpermutations)
    nullify(this%permsign)
    POP_SUB(permutations_nullify)
  end subroutine permutations_nullify

  subroutine permutations_copy (perm_in, perm_out)
    type(permutations_t), intent(inout) :: perm_in, perm_out

    PUSH_SUB(permutations_copy)

    perm_out%nn = perm_in%nn
    perm_out%npermutations = perm_in%npermutations
    perm_out%npairs = perm_in%npairs

    call loct_pointer_copy(perm_out%allpermutations,perm_in%allpermutations)
    call loct_pointer_copy(perm_out%permsign,perm_in%permsign)
    POP_SUB(permutations_copy)
  end subroutine permutations_copy

  subroutine permutations_end (this)
    type(permutations_t), intent(inout) :: this

    PUSH_SUB(permutations_end)
    SAFE_DEALLOCATE_P (this%allpermutations)
    SAFE_DEALLOCATE_P (this%permsign)
    POP_SUB(permutations_end)
  end subroutine permutations_end

end module permutations_m
