module liboct
  use global

  implicit none
  
  interface
    real(8) function oct_asinh(x)
      real(8), intent(in) :: x
    end function oct_asinh

    real(8) function oct_erf(x)
      real(8), intent(in) :: x
    end function oct_erf

    real(8) function oct_erfc(x)
      real(8), intent(in) :: x
    end function oct_erfc

    real(8) function oct_ylm(x, y, z, l, m)
      real(8), intent(in) :: x, y, z
      integer, intent(in) :: l, m
    end function oct_ylm

    integer function oct_parse_init(file_in, file_out)
      character(len=*), intent(in)  :: file_in, file_out
    end function oct_parse_init

    subroutine oct_parse_end()
    end subroutine oct_parse_end

    integer function oct_parse_isdef(name)
      character(len=*), intent(in) :: name
    end function oct_parse_isdef

    subroutine oct_parse_int(name, def, res)
      character(len=*), intent(in) :: name
      integer, intent(in)          :: def
      integer, intent(out)         :: res
    end subroutine oct_parse_int

    subroutine oct_parse_double(name, def, res)
      character(len=*), intent(in) :: name
      real(8), intent(in)          :: def
      real(8), intent(out)         :: res
    end subroutine oct_parse_double
    
    subroutine oct_parse_complex(name, def, res)
      character(len=*), intent(in) :: name
      complex(8), intent(in)       :: def
      complex(8), intent(out)      :: res
    end subroutine oct_parse_complex

    subroutine oct_parse_string(name, def, res)
      character(len=*), intent(in) :: name, def
      character(len=*), intent(out):: res
    end subroutine oct_parse_string

    integer function oct_parse_block_n(name)
      character(len=*), intent(in) :: name
    end function oct_parse_block_n

    subroutine oct_parse_block_int(name, l, c, res)
      character(len=*), intent(in) :: name
      integer, intent(in)          :: l, c
      integer, intent(out)         :: res
    end subroutine oct_parse_block_int

    subroutine oct_parse_block_double(name, l, c, res)
      character(len=*), intent(in) :: name
      integer, intent(in)          :: l, c
      real(8), intent(out)         :: res
    end subroutine oct_parse_block_double

    subroutine oct_parse_block_complex(name, l, c, res)
      character(len=*), intent(in) :: name
      integer, intent(in)          :: l, c
      complex(8), intent(out)      :: res
    end subroutine oct_parse_block_complex

    subroutine oct_parse_block_string(name, l, c, res)
      character(len=*), intent(in) :: name
      integer, intent(in)          :: l, c
      character(len=*), intent(out):: res
    end subroutine oct_parse_block_string

  end interface

  private :: oct_parse_string, oct_parse_block_string
contains
  
  ! logical is a FORTRAN type, so we emulate the routine with integers
  subroutine oct_parse_logical(name, def, res)
    character(len=*), intent(IN) :: name
    logical, intent(in) :: def
    logical, intent(out) :: res
    
    integer :: idef, ires
    
    idef = 0
    if(def) idef = 1
    
    call oct_parse_int(C_String(name), idef, ires)
    res = (ires .ne. 0)
    
  end subroutine oct_parse_logical

  ! to avoid errors
  subroutine oct_parse_str(name, def, res)
    character(len=*), intent(in)  :: name, def
    character(len=*), intent(out) :: res

    res = ""
    call oct_parse_string(C_string(name), C_string(def), res)
  end subroutine oct_parse_str

  subroutine oct_parse_block_logical(name, l, c, res)
    character(len=*), intent(IN) :: name
    integer, intent(in)          :: l, c
    logical, intent(out)         :: res
    
    integer :: ires
    
    call oct_parse_block_int(C_String(name), l, c, ires)
    res = (ires .ne. 0)
    
  end subroutine oct_parse_block_logical

  subroutine oct_parse_block_str(name, l, c, res)
    character(len=*), intent(in)  :: name
    integer, intent(in) :: l, c
    character(len=*), intent(out) :: res

    character(len=80) :: name1

    name1 = C_string(name)
    res = ""
    call oct_parse_block_string(name1, l, c, res)
  end subroutine oct_parse_block_str

  ! we need this routine to communicate strings to C
  character(len=80) function C_string(si) result(so)
    character(len=*), intent(IN) :: si
    
    integer :: i
    
    i = len(trim(si))
    so = ""
    so(1:i) = si(1:i)
    so(i+1:i+1) = achar(0)
  end function C_string

end module liboct
