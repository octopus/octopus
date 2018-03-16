#include "global.h"

module uuid_oct_m

  use global_oct_m
  use messages_oct_m
  use profiling_oct_m

  implicit none
  
  integer, parameter :: pd = 5
  integer, parameter :: sz = 8

  integer, parameter :: fp = selected_real_kind(p=pd+1)
  integer, parameter :: ip = selected_real_kind(r=pd)

  real(kind=fp),    parameter :: mltp = 10.0_fp**pd
  integer(kind=ip), parameter :: mask = 2**16-1

  private
  public :: &
    uuid_t

  public ::   &
    uuid_nil

  public ::        &
    assignment(=), &
    operator(==),  &
    operator(/=)
  
  public ::    &
    uuid_new,  &
    uuid_del,  &
    uuid_init, &
    uuid_hash, &
    uuid_char, &
    uuid_end
  
  type :: uuid_t
    private
    integer(kind=ip), dimension(sz) :: uuid = 0_ip
  end type uuid_t

  type(uuid_t), parameter :: uuid_nil = uuid_t()

  interface assignment(=)
    module procedure uuid_copy
  end interface assignment(=)

  interface operator(==)
    module procedure uuid_equl
  end interface operator(==)

  interface operator(/=)
    module procedure uuid_neql
  end interface operator(/=)

  interface uuid_new
    module procedure uuid_new_type
    module procedure uuid_new_copy
    module procedure uuid_new_char
  end interface uuid_new

  interface uuid_init
    module procedure uuid_init_type
    module procedure uuid_init_char
  end interface uuid_init

contains

  function uuid_new_type() result(this)
    type(uuid_t), pointer :: this

    PUSH_SUB(uuid_new_type)

    nullify(this)
    SAFE_ALLOCATE(this)
    call uuid_init(this)
    
    POP_SUB(uuid_new_type)
  end function uuid_new_type
  
  function uuid_new_copy(that) result(this)
    type(uuid_t), intent(in) :: that
    
    type(uuid_t), pointer :: this

    PUSH_SUB(uuid_new_copy)

    nullify(this)
    SAFE_ALLOCATE(this)
    call uuid_copy(this, that)
    
    POP_SUB(uuid_new_copy)
  end function uuid_new_copy
  
  function uuid_new_char(that) result(this)
    character(len=*), intent(in) :: that
    
    type(uuid_t), pointer :: this

    PUSH_SUB(uuid_new_char)

    nullify(this)
    SAFE_ALLOCATE(this)
    call uuid_init(this, that)
    
    POP_SUB(uuid_new_char)
  end function uuid_new_char
  
  subroutine uuid_del(this)
    type(uuid_t), pointer, intent(inout) :: this

    PUSH_SUB(uuid_del)

    if(associated(this))then
      call uuid_end(this)
      SAFE_DEALLOCATE_P(this)
    end if
    nullify(this)
    
    POP_SUB(uuid_del)
  end subroutine uuid_del
  
  subroutine uuid_init_type(this)
    type(uuid_t), intent(out) :: this

    real(kind=fp) :: flt
    integer       :: idx

    PUSH_SUB(uuid_init_type)

    do idx = 1, sz
      call random_number(flt)
      this%uuid(idx) = uuid_int16(flt)
    end do
    this%uuid(4) = ibclr(this%uuid(4), 12)
    this%uuid(4) = ibclr(this%uuid(4), 13)
    this%uuid(4) = ibset(this%uuid(4), 14)
    this%uuid(4) = ibclr(this%uuid(4), 15)
    call cpu_time(flt)
    this%uuid(5) = iand(mask, ieor(this%uuid(5), uuid_int16(flt)))
    this%uuid(5) = ibclr(this%uuid(5), 14)
    this%uuid(5) = ibset(this%uuid(5), 15)
        
    POP_SUB(uuid_init_type)
  end subroutine uuid_init_type
  
  subroutine uuid_init_char(this, that)
    type(uuid_t),     intent(out) :: this
    character(len=*), intent(in)  :: that

    PUSH_SUB(uuid_init_char)

    ASSERT(len_trim(that)>35)
    read(unit=that(1:4),   fmt="(z4.4)") this%uuid(1)
    read(unit=that(5:8),   fmt="(z4.4)") this%uuid(2)
    ASSERT(that(9:9)=='-')
    read(unit=that(10:13), fmt="(z4.4)") this%uuid(3)
    ASSERT(that(14:14)=='-')
    read(unit=that(15:18), fmt="(z4.4)") this%uuid(4)
    ASSERT(that(19:19)=='-')
    read(unit=that(20:23), fmt="(z4.4)") this%uuid(5)
    ASSERT(that(24:24)=='-')
    read(unit=that(25:28), fmt="(z4.4)") this%uuid(6)
    read(unit=that(29:32), fmt="(z4.4)") this%uuid(7)
    read(unit=that(33:36), fmt="(z4.4)") this%uuid(8)
        
    POP_SUB(uuid_init_char)
  end subroutine uuid_init_char
  
  elemental function uuid_int16(this) result(that)
    real(kind=fp), intent(in) :: this
    
    integer(kind=ip) :: that

    real(kind=fp) :: flt

    flt = 10.0_fp * fraction(this)
    flt = flt - real(floor(flt, kind=ip), kind=fp)
    flt =  real(floor(flt * mltp), kind=fp)
    that = iand(mask, int(flt))
        
  end function uuid_int16
  
  elemental function uuid_equl(this, that) result(equl)
    type(uuid_t), intent(in) :: this
    type(uuid_t), intent(in) :: that

    logical :: equl

    equl = all(iand(mask,this%uuid)==iand(mask, that%uuid))
        
  end function uuid_equl
  
  elemental function uuid_neql(this, that) result(neql)
    type(uuid_t), intent(in) :: this
    type(uuid_t), intent(in) :: that

    logical :: neql

    neql = any(iand(mask,this%uuid)/=iand(mask, that%uuid))
        
  end function uuid_neql
  
  elemental function uuid_hash(this, size) result(that)
    type(uuid_t),      intent(in) :: this
    integer, optional, intent(in) :: size

    integer :: that

    integer :: intl, inth

    inth = ieor(this%uuid(1), this%uuid(3))
    inth = ieor(inth, this%uuid(5))
    inth = ieor(inth, this%uuid(7))
    intl = ieor(this%uuid(2), this%uuid(4))
    intl = ieor(intl, this%uuid(6))
    intl = ieor(intl, this%uuid(8))
    that = ieor(ishft(inth, 15), intl)
    if(present(size)) that = modulo(that, size) + 1
        
  end function uuid_hash
  
  elemental function uuid_char(this) result(that)
    type(uuid_t), intent(in) :: this

    character(len=36) :: that

    write(unit=that, fmt="(2z4.4,'-',z4.4,'-',z4.4,'-',z4.4,'-',3z4.4)") this%uuid
        
  end function uuid_char
  
  elemental subroutine uuid_copy(this, that)
    type(uuid_t), intent(inout) :: this
    type(uuid_t), intent(in)    :: that

    this%uuid = iand(mask, that%uuid)
        
  end subroutine uuid_copy
  
  elemental subroutine uuid_end(this)
    type(uuid_t), intent(inout) :: this

    this%uuid = 0_ip
        
  end subroutine uuid_end
  
end module uuid_oct_m
