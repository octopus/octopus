!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This file contains modules io and fdf. The former takes care of the unit    !
!	management of any Fortran 90 code; the latter (which needs io) is     !
!	the fdf package.                                                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!










!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
! This module implements an interface to the FORTRAN logical unit             !
! system. Based on code by Richard Maine.                                     !
!                                                                             !
! Alberto Garcia, December 30, 1996                                           !
! Rewritten as a single subroutine                                            !
! with multiple entry points, March 7, 1998                                   !
!                                                                             !
! Converted into f90-module, A. Rubio 1999.                                   !
!                                                                             !
! Logical unit management. Units 0 to min_lun-1 are "reserved",               !
! since most of the "typical" files (output, etc) use them.                   !
!                                                                             !
! Logical units min_lun to min_max are managed by this module.                !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  module io

      
  implicit none




  integer, private, save  :: stdout=6, stderr=0
  integer, private, parameter :: min_lun=10, max_lun=99
  integer, private, parameter :: nunits=max_lun-min_lun+1
  integer, private   :: i, unit, lun, iostat

  logical, private, save  :: lun_is_free(min_lun:max_lun)=.true.
  logical, private   :: used, named, opened

  character(len=50), private :: filename
  character(len=11), private :: form









  contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Procedures                                                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Next subroutines find or change the standard units...                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine io_seterr(unit)

  integer, intent(inout) :: unit
  stderr = unit
  return

  end subroutine io_seterr




  subroutine io_setout(unit)

  integer, intent(inout) :: unit
  stdout = unit
  return

  end subroutine io_setout




  subroutine io_geterr(unit)

  integer, intent(inout) :: unit
  unit = stderr
  return

  end subroutine io_geterr




  subroutine io_getout(unit)

  integer, intent(inout) :: unit
  unit = stdout
  return

  end subroutine io_getout


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!











!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Logical unit management                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine io_assign(lun)


  integer, intent(INOUT) :: lun

! Looks for a free unit and assigns it to lun

  do lun= min_lun, max_lun
     if (lun_is_free(lun)) then
        inquire(unit=lun, opened=used, iostat=iostat)
        if (iostat .ne. 0) used = .true.
           lun_is_free(lun) = .false.
           if (.not. used) return
        endif
  enddo
  write(stderr,'(a)') 'No luns available in io_assign'
  stop 'LUN'

  end subroutine io_assign




  subroutine io_reserve(lun)


! Useful to specify that one needs to use a particular unit number
! For example, assume some legacy code expects to work with unit 15:
!
! call io_reserve(15)   ! this call at the beginning of the program
! ...
! open(15,....)

  integer, intent(INOUT) :: lun
  inquire(unit=lun, opened=used, iostat=iostat)
  if (iostat .ne. 0) used = .true.
  if (used) then
      write(stderr,'(a,i3,a)')                                     &
         'Cannot reserve unit',lun,'. Already connected'
      stop 'LUN'
  endif
      if (lun .ge. min_lun .and. lun .le. max_lun)                    &
                            lun_is_free(lun) = .false.

  return

  end subroutine io_reserve




  subroutine io_close(lun)


!     Use this routine instead of a simple close!!

  integer, intent(INOUT) :: lun
  close(lun)
  if (lun .ge. min_lun .and. lun .le. max_lun)                   &
                        lun_is_free(lun) = .true.
  return

  end subroutine io_close




  subroutine io_status


!     Prints a list of the connected logical units and the names of
!     the associated files


  write(stdout,'(a)') '******** io_status ********'
  do i = 0, max_lun
     inquire(i,opened=opened,named=named,name=filename,          &
             form=form,iostat=iostat)
     if (iostat .eq. 0) then
        if (opened) then
           if (named) then
               write(stdout,9000) i, form, filename
           else
               write(stdout,9000) i, form, 'No name available'
           endif
        endif
     else
        write(stdout,9000) i, 'Iostat error'
        endif
     enddo
     write(stdout,'(a)') '********           ********'

  9000 format(i4,5x,a,5x,a)

  return

  end subroutine io_status


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






  end module io


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
















!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! fdf module is a set of subroutines used to get data from an input file in   !
!	a comfortable way.                                                    !
!                                                                             !
! FDF (Flexible Data Format) routine package                                  !
!                                                                             !
! Jose Soler, Alberto Garcia Dec 1996                                         !
!                                                                             !
!                                                                             !
!                                                                             !
!     Notes:                                                                  !
!                                                                             !
!     This package implements fdf Version 0.6 (see file fdf.Standard)         !
!                                                                             !
!     User callable routines:                                                 !
!                                                                             !
!     fdf_init(filein,fileout)                                                !
!     fdf_shutdown                                                            !
!     fdf_enabled                                                             !
!     fdf_{integer,single,double,string,boolean} (label,default)              !
!     fdf_physical(label,default,unit)                                        !
!     fdf_block(label,io_unit)                                                !
!     fdf_defined(label)                                                      !
!                                                                             !
!     Implementation notes:                                                   !
!                                                                             !
!     This version needs the io module for logical unit management.           !
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  module fdf


  use io


  implicit none

!     I/O variables for the fdf package. 

  integer, private, parameter :: maxdepth = 5
  integer, private  :: ndepth, fdf_stack(maxdepth)

!     Unit numbers for input, output, error notification, and
!     debugging output (the latter active if fdf_debug is true)

  integer, private ::  fdf_in, fdf_out, fdf_err, fdf_log
  logical, private ::  fdf_debug, fdf_debug2, fdf_started

!     Line just read and parsing info

  character(len=132), private :: line
  integer, private, parameter :: maxntokens=50
  integer, private :: ntokens
  integer, private :: first(maxntokens), last(maxntokens)

  data fdf_started, fdf_debug, fdf_debug2 /.false.,.false.,.false./










  contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Procedures                                                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!









!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     New initialization for fdf. Simplified user interface using             !
!     the io package.                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine fdf_init(filein,fileout)


  implicit none

  character(len=*) filein, fileout

  integer debug_level

!     Prevent the user from opening two head files

  if (fdf_started) then 
      write(fdf_err,'(a)') 'FDF: Head file already set...'
      stop 'HEAD'
  endif

  call io_geterr(fdf_err)

  ndepth = 0

  call fdf_open(filein)
  call io_assign(fdf_out)
  open(unit=fdf_out,file=fileout,form='formatted',status='unknown')
  rewind(fdf_out)

  write(fdf_out,'(/,a,a,a,i3,/)')                                   &
        '#FDF: Opened ',filein, ' for input. Unit:',fdf_in

  fdf_started = .true.

  debug_level = fdf_integer('fdf-debug',0)
  call fdf_setdebug(debug_level)

  return

  end subroutine fdf_init


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Closes the 'head' file                                                  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine fdf_shutdown


  implicit none

  if (.not. fdf_started) return

  call fdf_refresh
  call io_close(fdf_in)
  fdf_started = .false.

  return

  end subroutine fdf_shutdown


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  logical function fdf_enabled()

  implicit none

  fdf_enabled = fdf_started
  return

  end function fdf_enabled


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine fdf_setdebug(level)

!     Debugging levels: 
!     level <=0: nothing
!     level  =1: standard
!     level >=2: exhaustive

   implicit none

   integer :: level

   if (level .le. 0) then

      if (fdf_debug) then
         call io_close(fdf_log)
         fdf_debug = .false.
      endif

   else

      if (.not. fdf_debug) then
         call io_assign(fdf_log)
         open(fdf_log,file='FDF.debug',form='formatted',status='unknown')
         rewind(fdf_log)
         fdf_debug = .true.
      endif

   endif
      
   fdf_debug2 = (level .ge. 2)

   return

   end subroutine fdf_setdebug


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine fdf_open(filename)

  implicit none

!     Opens a file for fdf processing.

  character(len=*) :: filename
  integer :: lun
  logical :: file_exists

  ndepth = ndepth + 1
  if (ndepth .gt. maxdepth) then
     write(fdf_err,'(a)') 'FDF: Too many nested fdf files...'
     stop 'DEPTH'
  endif

  if (leqi(filename,'stdin')) then
      lun = 5
      if (fdf_debug) write(fdf_log,'(a,i1,a)')                     &
              '--->Reading from Standard Input [depth:', ndepth,'] '

  else

      call io_assign(lun)
      inquire(file=filename,exist=file_exists)
      if (file_exists) then
      open(unit=lun,file=filename,status='old',form='formatted')
      rewind(lun)
      if (fdf_debug) write(fdf_log,'(a,i1,a,a50)')      &       
                 '--->Opened [depth:', ndepth,'] ', filename
      else
         write(fdf_err,'(a,a60)')                          &
              'FDF: Cannot open ',filename
      endif

  endif

      fdf_stack(ndepth) = lun
      fdf_in = lun
      
      return

  end subroutine fdf_open


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine fdf_close

  implicit none

!     Closes currently opened fdf file, except if it is the original one.

  if (ndepth .gt. 1) then
     call io_close(fdf_in)
     if (fdf_debug) write(fdf_log,'(a,i1,a)')'--->Closed [depth:', ndepth,']'

     ndepth = ndepth -1
     fdf_in = fdf_stack(ndepth)
  endif

  return

  end subroutine fdf_close


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine fdf_refresh
!
!     Closes all the open files in the stack (except the first).
!     Failure to do so would imply that the next Label is searched 
!     first in the 'deeper' files. fdf_locate calls fdf_refresh 
!     before doing anything. 

   implicit none
      
   integer i

   do i = ndepth, 1 , -1
     call fdf_close
  enddo
      
  return

  end subroutine fdf_refresh


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine fdf_parse


!     Processes the input line looking for meaningful tokens.

  implicit none

  logical intoken, instring

  integer c
  integer stringdel

!     Character statement functions

  integer i
  logical isdigit, isupper, islower, isalpha, isalnum, isextra, istokch
  logical iscomment, isdelstr, isspecial

      isdigit(i) = (i .ge. 48) .and. (i .le. 57)
      isupper(i) = (i .ge. 65) .and. (i .le. 90)
      islower(i) = (i .ge. 97) .and. (i .le. 122)
      isalpha(i) = isupper(i) .or. islower(i)
      isalnum(i) = isdigit(i) .or. isalpha(i)

!     Extra characters allowed in tokens:  $ % * + & - . / @ ^ _ ~
      isextra(i) = ((i .ge. 36) .and. (i .le. 38))                     &
                 .or. (i .eq. 42) .or. (i .eq. 43)                     &
                 .or. (i .eq. 45) .or. (i .eq. 46)                     &
                 .or. (i .eq. 47) .or. (i .eq. 64) .or. (i .eq. 94)    &
                 .or. (i .eq. 95) .or. (i .eq. 126)

      istokch(i) = isalnum(i) .or. isextra(i)

!     Comments are signaled by:  !  #  ; 
      iscomment(i) = (i.eq.33) .or. (i.eq.35) .or. (i.eq.59)

!     String delimiters:
      isdelstr(i) = (i.eq.34) .or. (i.eq.39) .or. (i.eq.96)

!     Special characters which are tokens by themselves: <
      isspecial(i) = (i.eq.60)



      intoken = .false.
      instring = .false.
      ntokens = 0
      stringdel = 0
      
      do i = 1, len(line)
         c = ichar(line(i:i))

         if (iscomment(c)) then
! possible comment...
            if (instring) then
               last(ntokens) = i
            else
               goto 1000
            endif

         else if (istokch(c)) then
! character allowed in a token...
            if (.not. intoken) then
               intoken = .true.
               ntokens = ntokens+1
               first(ntokens) = i
            endif
            last(ntokens) = i

         else if (isspecial(c)) then
! character that forms a token by itself...
            if (.not. instring) then
               ntokens=ntokens+1
               first(ntokens) = i
               intoken = .false.
            endif
            last(ntokens) = i

         else if (isdelstr(c)) then
! string delimiter... make sure it is the right one before
! closing the string.
! If we are currently in a token, the delimiter is appended to it.

            if (instring) then
               if (c.eq.stringdel) then
                  instring = .false.
                  intoken = .false.
                  stringdel = 0
               else
                  last(ntokens) = i
               endif
            else
               if (intoken) then
                  last(ntokens) = i
               else
                  instring = .true.
                  stringdel = c
                  intoken = .true.
                  ntokens = ntokens+1
                  first(ntokens) = i+1
                  last(ntokens) = i+1
               endif
            endif

         else
! token delimiter...

            if (instring) then
               last(ntokens) = i
            else
               if (intoken) intoken=.false.
            endif
         endif

      enddo

 1000 continue

      if (fdf_debug2) then
         write(fdf_log,*) '            ',  ntokens, ' tokens:'
         do i=1,ntokens
            write(fdf_log,*) '                 ',              &
                             '|',line(first(i):last(i)),'|'
         enddo
      endif


  return

  end subroutine fdf_parse


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  integer function fdf_search(label)

!     Performs a case-and-punctuation-insensitive search for 'label'
!     among the tokens in a line.

  implicit none

  character(len=*) label
      integer i

      fdf_search = 0
      do i = 1, ntokens
         if (labeleq(label,line(first(i):last(i)))) then
            fdf_search = i
            return
         endif
      enddo

      return
  end function fdf_search


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  character*(80) function fdf_string(label,default)

!     Returns a string associated with label label, or default if label
!     is not found in the fdf file.

      implicit none

      character*(*) label
      character*(*) default
!
      fdf_string = default

      if (.not. fdf_locate(label)) then
         write(fdf_out,'(a,5x,a,5x,a)') label, default,'# Default value'
         return
      endif
!
!     From the second token up...
!
      fdf_string = line(first(2):last(ntokens))
      write(fdf_out,'(a,5x,a)') label, fdf_string

      return

  end function fdf_string


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      logical function fdf_boolean(label,default)
!
!     Returns true if label appears by itself or in the form
!     label {Yes,true,.true.,T} (case insensitive).
!
!     Returns false if label appears in the form
!     label {No,false,.false.,F} (case insensitive).
!
!     If label is not found in the fdf file, fdf_boolean returns the 
!     logical variable default.
!
!     implicit none

      character*(*) label
      logical default

      character valstr*40
!
      fdf_boolean = default

      if (.not. fdf_locate(label)) then
         write(fdf_out,'(a,5x,l10,5x,a)')                      &
                        label, default, '# Default value'
         return
      endif

!
!     If the label appears by itself, we interpret it as .true.
!
      if (ntokens .eq. 1) then
         fdf_boolean = .true.
         write(fdf_out,'(a,5x,l10,5x,a)') label, fdf_boolean,          &
                                          '# Label by itself'
         return
      endif
!
!     Look for second word
!
      valstr=line(first(2):last(2))
!
      if (leqi(valstr,'yes') .or. leqi(valstr,'true') .or.              &
          leqi(valstr,'.true.') .or. leqi(valstr,'t') .or.              &
          leqi(valstr,'y'))       then
         
         fdf_boolean = .true.
         write(fdf_out,'(a,5x,l10)') label, fdf_boolean

      else if (leqi(valstr,'no') .or. leqi(valstr,'false') .or.         &
               leqi(valstr,'.false.') .or. leqi(valstr,'f') .or.        &
               leqi(valstr,'n'))       then

         fdf_boolean = .false.
         write(fdf_out,'(a,5x,l10)') label, fdf_boolean

      else

         write(fdf_err,*)                                               &
              'FDF_BOOLEAN: Unexpected fdf logical value ',             &
              label, ' = ', valstr
         stop 

      endif

  return

  end function fdf_boolean


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  integer function fdf_integer(label,default)


!     Returns an integer associated with label, or default if label
!     is not found in the fdf file.

      implicit none

      character*(*) label
      integer default

      character*10 fmtstr

      fdf_integer = default

      if (.not. fdf_locate(label)) then
         write(fdf_out,'(a,5x,i10,5x,a)')                            &
                         label, default, '# Default value'
         return
      endif

      if (ntokens.eq.1) then
         write(fdf_err,*) 'FDF_INTEGER: No value for ', label
         stop
      endif

      write(fmtstr,9000) last(2)-first(2)+1
 9000 format('(i',i2.2,')')
      read(line(first(2):last(2)),fmt=fmtstr) fdf_integer
      write(fdf_out,'(a,5x,i20)') label, fdf_integer

  return

  end function fdf_integer


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  real function fdf_single(label,default)

!     Returns a single precision value associated with label label, 
!     or default if label is not found in the fdf file.

      implicit none

      character*(*) label
      real default
!
      character*10 fmtstr
!
      fdf_single = default

      if (.not. fdf_locate(label)) then
         write(fdf_out,'(a,5x,g20.10,5x,a)')                            &
                         label, default, '# Default value'
         return
      endif

      if (ntokens.eq.1) then
         write(fdf_err,*) 'FDF_SINGLE: No value for ', label
         stop
      endif
      write(fmtstr,9000) last(2)-first(2)+1
 9000 format('(g',i2.2,'.0)')
      read(line(first(2):last(2)),fmt=fmtstr) fdf_single
      write(fdf_out,'(a,5x,g20.10)') label, fdf_single


  return

  end function fdf_single


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  double precision function fdf_double(label,default)

!     Returns a double precision value associated with label label, 
!     or default if label is not found in the fdf file.

  implicit none

      character*(*) label
      double precision default
!
      character*10 fmtstr
!
      fdf_double = default

      if (.not. fdf_locate(label)) then
         write(fdf_out,'(a,5x,g20.10,5x,a)')                        &
                         label, default, '# Default value'
         return
      endif

      if (ntokens.eq.1) then
         write(fdf_err,*) 'FDF_DOUBLE: No value for ', label
         stop
      endif
      write(fmtstr,9000) last(2)-first(2)+1
 9000 format('(g',i2.2,'.0)')
      read(line(first(2):last(2)),fmt=fmtstr) fdf_double
      write(fdf_out,'(a,5x,g20.10)') label, fdf_double


  return

  end function fdf_double


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  double precision function fdf_physical(label,default,defunit)

!     Returns a double precision value associated with label label, 
!     or default if label is not found in the fdf file. Converts
!     the units to defunit.

  implicit none

      character*(*) label, defunit
      double precision default

      character unitstr*10
      double precision value

      character*10 fmtstr
!
      fdf_physical = default

      if (.not. fdf_locate(label)) then
         write(fdf_out,'(a,5x,g20.10,1x,a,5x,a)')                     &
                         label, default, defunit, '# Default value'
         return
      endif

      if (ntokens.eq.1) then
         write(fdf_err,*) 'FDF_PHYSICAL: No value for ', label
         stop
      endif
      write(fmtstr,9000) last(2)-first(2)+1
 9000 format('(g',i2.2,'.0)')
      read(line(first(2):last(2)),fmt=fmtstr) value
      fdf_physical = value
!
!     Look for unit
!
      if (ntokens.eq.2) then
         write(fdf_err,*) 'FDF_PHYSICAL: No unit specified for ', label
         stop
      endif

      unitstr=line(first(3):last(3))
      if (.not. leqi(unitstr,defunit))                                   &
           fdf_physical = value * fdf_convfac(unitstr,defunit)
      write(fdf_out,'(a,5x,g20.10,1x,a10)')                              &
                     label, fdf_physical, defunit
      write(fdf_out,'(a,a,5x,g20.10,1x,a10)')                            &
                     '# Above item originally: ', label, value, unitstr

  return

  end function fdf_physical


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  logical function fdf_block(label,unit)

!     Returns "true" and the unit number of the file from which to read
!     the contents of a block if "label" is associated with a block, and
!     false if not (unit is set to -1 in this case).

      implicit none
      character*(*) label
      integer unit
      
      character*50 token1, filename
      integer iless

!      logical fdf_locate, leqi
!      integer fdf_search

      fdf_block = .false.
      unit = -1
      if (.not. fdf_locate(label)) return

      token1 = line(first(1):last(1))
      if (.not. leqi(token1,'%block')) then
         write(fdf_err,*) 'FDF_BLOCK: Not a block:',label
!
!        Return instead of stopping
!
         return
      endif

      iless = fdf_search('<')
      if ((iless .ne. 0) .and. (ntokens .gt. iless)) then
!
!           Read block from file
!
         filename = line(first(iless+1):last(iless+1))
         if (fdf_debug) write(fdf_log,'(2a)')                       &
              '*Reading block from file', filename
         call fdf_open(filename)
         if (fdf_search('%dump') .ne. 0) call fdf_dumpfile(label)
         fdf_block = .true.
         unit = fdf_in
         return
      endif
!
!     Standard block in fdf file. Dump contents
!
      call fdf_dumpblock(label)
      fdf_block = .true.
      unit = fdf_in
      
  return

  end function fdf_block


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  logical function fdf_locate(label)


!     Searches for label in the fdf hierarchy. If it appears and it
!     is not part of a comment, the function returns .true. and leaves
!     the file positioned at the next line. Otherwise, it returns .false.
!
!     It supports two kinds of "include" files:
!
!     %include filename  
!     Indicates an unconditional opening of filename for 
!     further fdf processing.
!
!     Label1 Label2 ... < filename  
!     Indicates that filename should be opened only when 
!     searching for any of the labels indicated.
!     'filename' should be an fdf file.
!
      implicit none
      
      character*(*) label

      character*60 token1, filename
      integer ilabel, iless
!
      call fdf_refresh
      if (fdf_debug) write(fdf_log,'(/,a,1x,a)') 'Looking for ', label

      fdf_locate = .false.

      rewind(fdf_in)

 10   continue

      if (.not. fdf_getline()) then
         if (ndepth .gt. 1) then
            call fdf_close
            goto 10
         endif
         if (fdf_debug) write(fdf_log,'(a,1x,a)') '*Did not find ', label
         return
      endif
!
      if (ntokens .eq. 0) goto 10
!
      token1 = line(first(1):last(1))
!
      if (leqi(token1,'%include')) then
!
!        Include file
!
         if (ntokens .eq. 1) then
            write(fdf_err,*) 'FDF: No valid filename after %include'
            stop
         endif
         filename = line(first(2):last(2))
         call fdf_open(filename)
         goto 10
      endif

      ilabel = fdf_search(label)
      if (ilabel .ne. 0) then
!
!        Label found...
!
         if (leqi(token1,'%block')) then
            fdf_locate = .true.
            if (fdf_debug) write(fdf_log,'(a,1x,a)')'*Found ', label
            return
         endif

         iless = fdf_search('<')
         if ((iless .ne. 0) .and. (ntokens .gt. iless)) then
!
!           Continue search in other file
!
            filename = line(first(iless+1):last(iless+1))
            call fdf_open(filename)
            goto 10
         endif
!
!        If we reach this point we must be dealing with a line
!        of the form Label Value. But we are not interested if
!        the string appears in the Value section
!
         if (ilabel .eq. 1) then
            fdf_locate = .true.
            if (fdf_debug) write(fdf_log,'(a,1x,a)') '*Found ', label
            return
         else
            goto 10
         endif

      else

         goto 10

      endif


  end function fdf_locate


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  logical function fdf_defined(label)

  implicit none
      
  character(len=*) label
      
  fdf_defined = fdf_locate(label)

  return

  end function fdf_defined


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  logical function fdf_getline()

  implicit none

  read(fdf_in,end=100,err=100,fmt='(a)') line
      fdf_getline = .true.
      if (fdf_debug2) write(fdf_log,'(a,a76)') '> ', line
      call fdf_parse
      return

 100  continue
      fdf_getline = .false.

  return

  end function fdf_getline


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  double precision function fdf_convfac( from, to )

! Returns conversion factor between a subset of physical units
! Written by J.M.Soler. Dec 96.
! Modified by Alberto Garcia, Jan 97.

      IMPLICIT      NONE
      CHARACTER*(*) FROM, TO
      INTEGER       IU, IFROM, ITO, NU

      PARAMETER ( NU = 54 )
      CHARACTER         DIM(NU)*10, NAME(NU)*10
      DOUBLE PRECISION  UNIT(NU)


!     We allow case variations in the units. This could be dangerous
!     (meV --> MeV !!) in real life, but not in this restricted 
!     field.
     


      DATA (DIM(IU), NAME(IU), UNIT(IU), IU=1,10) /          &
        'mass  ', 'Kg      ', 1.D0,                          &
        'mass  ', 'g       ', 1.D-3,                         &
        'mass  ', 'amu     ', 1.66054D-27,                   &
        'length', 'm       ', 1.D0,                          &
        'length', 'nm      ', 1.D-9,                         &
        'length', 'Ang     ', 1.D-10,                        &
        'length', 'Bohr    ', 0.529177D-10,                  &
        'time  ', 's       ', 1.D0,                          &
        'time  ', 'fs      ', 1.D-15,                        &
        'energy', 'J       ', 1.D0/
      DATA (DIM(IU), NAME(IU), UNIT(IU), IU=11,20) /          &
        'energy', 'erg     ', 1.D-7,                          &
        'energy', 'eV      ', 1.60219D-19,                    &
        'energy', 'meV     ', 1.60219D-22,                    &
        'energy', 'Ry      ', 2.17991D-18,                    &
        'energy', 'mRy     ', 2.17991D-21,                    &
        'energy', 'Hartree ', 4.35982D-18,                    &
        'energy', 'K       ', 1.38066D-23,                    &
        'energy', 'kcal/mol', 6.94780D-21,                    &
        'force ', 'N       ', 1.D0,                           &
        'force ', 'eV/Ang  ', 1.60219D-9/
      DATA (DIM(IU), NAME(IU), UNIT(IU), IU=21,30) /          &
        'force ', 'Ry/Bohr ', 4.11943D-8,                     &
        'length  ', 'cm      ', 1.d-2,                        &
        'time    ', 'ps      ', 1.d-12,                       &
        'time    ', 'ns      ', 1.d-9,                        &
        'energy  ', 'mHartree', 4.35982D-21,                  &
        'energy  ', 'kJ/mol  ', 1.6606d-21,                   &
        'energy  ', 'Hz      ', 6.6262d-34,                   &
        'energy  ', 'THz     ', 6.6262d-22,                   &
        'energy  ', 'cm-1    ', 1.986d-23,                    &
        'energy  ', 'cm^-1   ', 1.986d-23/
      DATA (DIM(IU), NAME(IU), UNIT(IU), IU=31,40) /          &
        'pressure', 'Pa      ', 1.d0,                         &
        'pressure', 'MPa     ', 1.d6,                         &
        'pressure', 'GPa     ', 1.d9,                         &
        'pressure', 'atm     ', 1.01325d5,                    &
        'pressure', 'bar     ', 1.d5,                         &
        'pressure', 'Mbar    ', 1.d11,                        &
        'charge  ', 'C       ', 1.d0,                         &
        'charge  ', 'e       ', 1.602177d-19,                 &
        'dipole  ', 'C*m     ', 1.d0,                         &
        'dipole  ', 'D       ', 3.33564d-30/
      DATA (DIM(IU), NAME(IU), UNIT(IU), IU=41,50) /          &
        'dipole  ', 'debye   ', 3.33564d-30,                  &
        'dipole  ', 'e*Bohr  ', 8.47835d-30,                  &
        'dipole  ', 'e*Ang   ', 1.602177d-29,                 &
        'energy  ', 'cm**-1    ', 1.986d-23,                  &
        'pressure', 'Ry/Bohr**3', 1.47108d13,                 &
        'pressure', 'eV/Ang**3 ', 1.60219d11,                 &
        'MomInert', 'Kg*m**2   ', 1.d0,                       &
        'MomInert', 'Ry*fs**2  ', 2.17991d-48,                &
        'Efield  ', 'V/m       ', 1.d0,                       &
        'Efield  ', 'V/nm      ', 1.d9 /
      DATA (DIM(IU), NAME(IU), UNIT(IU), IU=51,54) /          &
        'Efield  ', 'V/Ang     ', 1.d10,                      &
        'Efield  ', 'V/Bohr    ', 1.8897268d10,               &
        'Efield  ', 'Ry/Bohr/e ', 2.5711273d11,               &
        'Efield  ', 'Har/Bohr/e', 5.1422546d11 /
!
      IFROM = 0
      ITO   = 0
      DO 10 IU = 1,NU
        IF (leqi(NAME(IU),FROM)) IFROM = IU
        IF (leqi(NAME(IU),TO))   ITO   = IU
   10 CONTINUE

      IF (IFROM .EQ. 0) THEN
        WRITE(6,*) 'FDF_CONVFAC: Unknown unit = ', FROM
        STOP
      ENDIF
      IF (ITO .EQ. 0) THEN
        WRITE(6,*) 'FDF_CONVFAC: Unknown unit = ', TO
        STOP
      ENDIF

      IF (leqi(DIM(IFROM),DIM(ITO))) THEN
        FDF_CONVFAC = UNIT(IFROM) / UNIT(ITO)
      ELSE
        WRITE(6,*)                                                      &
          'FDF_CONVFAC: Unit''s physical dimensions don''t match: ',    &
           FROM, ', ', TO
        STOP
      ENDIF

  end function fdf_convfac


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine fdf_dumpblock(label)
     
!     Dumps block contents starting at the current line
     
      implicit none

      character*(*) label

      integer i, lblock
      character*60 token1

!      logical fdf_getline, leqi
!      integer fdf_search

      write(fdf_out,'(/,a79)') line
      lblock = 0
 120  continue
      if (fdf_getline()) then
         lblock = lblock + 1
         write(fdf_out,'(a79)') line
         token1 = line(first(1):last(1))
         if (.not. leqi(token1,'%endblock')) goto 120
      else
         write(fdf_err,'(a,a,a)')'FDF_LOCATE: Block ', label, ' does not end!'
         stop 'FDF'
      endif
      write(fdf_out,*)
!     
!     Sanity check (optional construct %endblock [ Label [Label] ])
!     
      if ((ntokens .gt. 1) .and. (fdf_search(label) .eq. 0)) then
         write(fdf_err,'(a,a,a)')'FDF_LOCATE: Block ', label, ' does not end!'
         stop 'FDF'
      endif
!
!     Backspace the lines read
!
      do i=1,lblock
         backspace(fdf_in)
      enddo

  return

  end subroutine fdf_dumpblock


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     
  subroutine fdf_dumpfile(label)
     
!     Dumps the contents of a file to fdf_out.
!     The lines are embedded in a %block ... %endblock pair.
!     
      implicit none

      character*(*) label
      character form*30
      integer length
!     
!     Build the right format
!     
      call chrlen(label,0,length)
      write(form,'(a,i2.2,a)') '(a,a',length,',10x,a)'

      write(fdf_out,*)
      write(fdf_out,form) '%block ', label,'# Originally in include file' 
!     
      rewind(fdf_in)
 10   continue
      if (fdf_getline()) then
         write(fdf_out,'(a79)') line
         goto 10
      endif

      write(fdf_out,form) '%endblock ', label, '# Originally in include file' 
      write(fdf_out,*)
      rewind(fdf_in)
!     
  return

  end subroutine fdf_dumpfile


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


   logical function labeleq(s1,s2)

!     Compares s1 and s2 without regard for case, or appearance
!     of '_', '.', '-'.

      implicit none

      character*(*) s1, s2
      character*80 n1, n2
!      logical leqi

      call fdf_pack(s1,n1)
      call fdf_pack(s2,n2)
      labeleq=leqi(n1,n2)
      if (fdf_debug) then
         if (labeleq .and. .not. leqi(s1,s2))                  &
              write(fdf_log,'(a,/,a,/,a)')                     &
              '--------- Considered equivalent:', s1, s2
      endif

  return

  end function labeleq


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine fdf_pack(s,n)

  implicit none
      character*(*) s, n
!
!     Removes occurrences of '_ .-'  from s1
!
      character*1 c
      integer i, j
      logical issep
      issep(i) = (i.eq.95) .or. (i.eq.46) .or. (i.eq.45)

      n = ' '
      j = 0
      do i = 1, len(s)
         c = s(i:i)
         if (.not.issep(ichar(c))) then
            j = j+1
            n(j:j) = c
         endif
      enddo

  return

  end subroutine fdf_pack


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine chrlen(string,nchar,lchar)


!  CHRLEN accepts a STRING of NCHAR characters and returns LCHAR,
!  the length of the string up to the last nonblank, nonnull.
     
      implicit none

      CHARACTER CHAR*1
      CHARACTER STRING*(*)
      integer nchar,lchar
!
      integer ncopy, i
!
      NCOPY=NCHAR
      IF(NCOPY.LE.0)NCOPY=LEN(STRING)
!
      DO 10 I=1,NCOPY
        LCHAR=NCOPY+1-I
        IF(STRING(LCHAR:LCHAR).NE.' '.AND.STRING(LCHAR:LCHAR).NE.CHAR(0))RETURN
 10     CONTINUE
      LCHAR=0

  return

  end subroutine chrlen


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  subroutine chrcap(STRING,NCHAR)
!
!  CHRCAP accepts a STRING of NCHAR characters and replaces
!  any lowercase letters by uppercase ones.
!
      implicit none

      CHARACTER CHAR*1
      integer nchar, ncopy, i, itemp
      LOGICAL   LGE
      LOGICAL   LLE
      CHARACTER STRING*(*)
!
      NCOPY=NCHAR
      IF(NCOPY.LE.0)NCOPY=LEN(STRING)
      DO 10 I=1,NCOPY
!
        IF(LGE(STRING(I:I),'a').AND.LLE(STRING(I:I),'z'))THEN
          ITEMP=ICHAR(STRING(I:I))+ICHAR('A')-ICHAR('a')
          STRING(I:I)=CHAR(ITEMP)
          ENDIF
10      CONTINUE

  return 

  end subroutine chrcap


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  logical function leqi(strng1,strng2)


!  Case-insensitive lexical equal-to comparison

  implicit none

  character :: s1,s2 
  character(len=*) :: strng1,strng2

  integer :: len1, len2, lenc, i

      LEN1=LEN(STRNG1)
      LEN2=LEN(STRNG2)
      LENC=MIN(LEN1,LEN2)

      LEQI=.FALSE.
      DO 10 I=1,LENC
        S1=STRNG1(I:I)
        S2=STRNG2(I:I)
        CALL CHRCAP(S1,1)
        CALL CHRCAP(S2,1)
        IF(S1.NE.S2)RETURN
10      CONTINUE
 
      IF(LEN1.GT.LENC.AND.STRNG1(LENC+1:LEN1).NE.' ')RETURN
      IF(LEN2.GT.LENC.AND.STRNG2(LENC+1:LEN2).NE.' ')RETURN
      LEQI=.TRUE.

  return 


  end function leqi









  end module fdf


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
