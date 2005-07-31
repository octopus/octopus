!! Copyright (C) 2002 M. Marques, A. Castro, A. Rubio, G. Bertsch
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

#include "global.h"

module hgh
  ! For information about the Hartwinger-Goedecker-Hutter pseudopotentials, take a look at:
  !  (1) S. Goedecker, M. Teter and J. Hutter, Phys. Rev. B 54, 1703 (1996).
  !  (2) C. Hartwinger, S. Goedecker and J. Hutter, Phys. Rev. B 58, 3641 (1998).
  use global
  use messages
  use io
  use atomic
  use logrid
  use lib_oct

  implicit none

  private
  public  :: hgh_type, hgh_init, hgh_process, hgh_debug, hgh_end

  ! Next data type contains:
  !   (a) the pseudopotential parameters, as read from a *.hgh file,
  !   (b) auxiliary intermediate functions, to store stuff before passing it to the "ps" variable.
  type hgh_type
     ! HGH parameters.
     character(len=5) :: atom_name
     integer          :: z_val
     FLOAT         :: rlocal
     FLOAT         :: rc(0:3)
     FLOAT         :: c(1:4)
     FLOAT         :: h(0:3, 1:3, 1:3)
     FLOAT         :: k(0:3, 1:3, 1:3)

     type(valconf)    :: conf
     integer          :: l_max     ! Maximum l for the Kleinmann-Bylander component.

     ! Logarithmic grid parameters
     type(logrid_type) :: g

     FLOAT, pointer :: vlocal(:) ! Local potential
     FLOAT, pointer :: kb(:,:,:) ! KB projectors
     FLOAT, pointer :: kbr(:)    ! KB radii
     FLOAT, pointer :: rphi(:,:), eigen(:)
  end type hgh_type

  FLOAT, parameter :: eps = CNST(1.0e-8)

  interface vlocalr
     module procedure vlocalr_scalar, vlocalr_vector
  end interface
  interface projectorr
     module procedure projectorr_scalar, projectorr_vector
  end interface

contains

  subroutine hgh_init(psp, filename)
    type(hgh_type), intent(inout)     :: psp
    character(len=*), intent(in)      :: filename

    integer :: iunit, i
    logical :: found
    character(len=256) :: filename2

    call push_sub('hgh.hgh_init')

    filename2 =  filename // '.hgh'
    inquire(file=filename2, exist=found)
    if(.not.found) then
       filename2 = trim(conf%share) // "/PP/HGH/" // filename // ".hgh"
       inquire(file=filename2, exist=found)
       if(.not.found) then
          message(1) = "Pseudopotential file '" // trim(filename) // ".hgh' not found!"
          call write_fatal(1)
       end if
    end if

    message(1) = "Info: Reading pseudopotential from file:"
    write(message(2), '(6x,3a)') "'", trim(filename2), "'"
    call write_info(2)

    iunit = io_open(filename2, action='read', form='formatted', status='old')
    i = load_params(iunit, psp)
    if(i .ne. 0) then
       message(1) = 'Error reading hgh file!'
       call write_fatal(1)
    endif
    call io_close(iunit)

    ! Finds out psp%l_max. The most special cases are H, He, Li_sc and Be_sc, where psp%l_max = -1.
    psp%l_max = 0
    do while(psp%rc(psp%l_max)>CNST(0.01))
       psp%l_max = psp%l_max + 1
    end do
    psp%l_max = psp%l_max - 1

    ! Initializes the logarithmic grid. Parameters are hard-coded.
    call logrid_init(psp%g, CNST(3.0e-2), CNST(4.0e-4), 431)

    ! Allocation of stuff.
    allocate(psp%vlocal(1:psp%g%nrval))
    psp%vlocal = M_ZERO
    if(psp%l_max >= 0) then
       allocate(psp%kbr(0:psp%l_max), &
            psp%kb(1:psp%g%nrval, 0:psp%l_max, 1:3))
       psp%kbr = M_ZERO; psp%kb = M_ZERO
    endif
    allocate(psp%rphi (psp%g%nrval, psp%conf%p), &
         psp%eigen(psp%conf%p))
    psp%rphi = M_ZERO; psp%eigen = M_ZERO

    call pop_sub()
  end subroutine hgh_init

  subroutine hgh_end(psp)
    type(hgh_type), intent(inout) :: psp

    call push_sub('hgh.hgh_end')

    deallocate(psp%vlocal, psp%rphi, psp%eigen)
    if(associated(psp%kb)) deallocate(psp%kb, psp%kbr)
    call logrid_end(psp%g)

    call pop_sub()
  end subroutine hgh_end

  subroutine hgh_process(psp)
    type(hgh_type), intent(inout) :: psp

    integer :: l, i

    call push_sub('hgh.hgh_process')

    ! Fixes the local potential
    psp%vlocal(1:psp%g%nrval) = vlocalr(psp%g%rofi, psp)

    ! And the projectors
    do l = 0, psp%l_max
       do i = 1, 3
          psp%kb(1:psp%g%nrval, l, i) = projectorr(psp%g%rofi, psp, i, l)
       enddo
    enddo

    ! get the pseudoatomic eigenfunctions (WARNING: This is not correctly done yet: "some" wavefunctions
    ! are obtained, but not the real ones!!!
    call solve_schroedinger(psp)

    ! Define the KB-projector cut-off radii
    call get_cutoff_radii(psp)

    call pop_sub()
  end subroutine hgh_process

  function load_params(unit, params)
    integer, intent(in)             :: unit        ! where to read from
    type(hgh_type), intent(out)     :: params      ! obvious
    integer                         :: load_params ! 0 if success,
    ! 1 otherwise.

    integer :: i, iostat, j, k
    character(len=VALCONF_STRING_LENGTH) :: line

    call push_sub('hgh.load_params')

    ! Set initially everything to zero.
    params%c(1:4) = M_ZERO; params%rlocal = M_ZERO;
    params%rc = M_ZERO; params%h = M_ZERO; params%k = M_ZERO

    ! get valence configuration
    read(unit,'(a)') line
    call read_valconf(line, params%conf)

    ! Reads the file in a hopefully smart way
    iostat = 1; j = 5
    read(unit,'(a)') line
    do while((iostat .ne. 0) .and. (j > 0))
       j = j - 1
       read(line, *, iostat=iostat) params%atom_name, params%z_val, params%rlocal, params%c(1:j)
    end do
    if(j<1) read(line, *, iostat=iostat) params%atom_name, params%z_val, params%rlocal
    if( iostat.ne.0 ) then
       load_params = 1
       call pop_sub(); return
    endif

    read(unit,'(a)', iostat = iostat) line
    if(iostat .ne. 0) then
       load_params = 0
       call pop_sub(); return
    endif
    iostat = 1; j = 4
    do while((iostat .ne. 0) .and. (j > 0))
       j = j - 1
       read(line, *, iostat=iostat) params%rc(0), (params%h(0, i, i), i = 1, j)
    enddo
    if(j < 0) then
       load_params = 2
       call pop_sub(); return
    endif

    kloop: do k = 1, 3
       read(unit, '(a)', iostat = iostat) line
       if(iostat .ne. 0) exit kloop
       iostat = 1; j = 4
       do while((iostat .ne. 0) .and. (j > 0))
          j = j - 1
          read(line, *, iostat = iostat) params%rc(k), (params%h(k, i, i), i = 1, j)
       enddo
       if(params%rc(k) == M_ZERO) exit kloop
       read(unit, '(a)') line
       iostat = 1; j = 4
       do while((iostat .ne. 0) .and. (j>0))
          j = j - 1
          read(line, *, iostat = iostat) (params%k(k, i, i), i = 1, 3)
       enddo
    enddo kloop

    ! Fill in the rest of the parameters matrices...
    ! Fill in the rest of the parameters matrices...
    params%h(0, 1, 2) = -M_HALF      * sqrt(M_THREE/M_FIVE)            * params%h(0, 2, 2)
    params%h(0, 1, 3) =  M_HALF      * sqrt(M_FIVE/CNST(21.0))         * params%h(0, 3, 3)
    params%h(0, 2, 3) = -M_HALF      * sqrt(CNST(100.0)/CNST(63.0))    * params%h(0, 3, 3)
    params%h(1, 1, 2) = -M_HALF      * sqrt(M_FIVE/M_SEVEN)            * params%h(1, 2, 2)
    params%h(1, 1, 3) =  M_ONE/M_SIX * sqrt(CNST(35.0)/CNST(11.0))     * params%h(1, 3, 3)
    params%h(1, 2, 3) = -M_ONE/M_SIX * (CNST(14.0) / sqrt(CNST(11.0))) * params%h(1, 3, 3)
    params%h(2, 1, 2) = -M_HALF      * sqrt(M_SEVEN/M_NINE)            * params%h(2, 2, 2)
    params%h(2, 1, 3) =  M_HALF      * sqrt(CNST(63.0)/CNST(143.0))    * params%h(2, 3, 3)
    params%h(2, 2, 3) = -M_HALF      * (CNST(18.0)/sqrt(CNST(143.0)))  * params%h(2, 3, 3)

    params%k(0, 1, 2) = -M_HALF      * sqrt(M_THREE/M_FIVE)            * params%k(0, 2, 2)
    params%k(0, 1, 3) =  M_HALF      * sqrt(M_FIVE/CNST(21.0))         * params%k(0, 3, 3)
    params%k(0, 2, 3) = -M_HALF      * sqrt(CNST(100.0)/CNST(63.0))    * params%k(0, 3, 3)
    params%k(1, 1, 2) = -M_HALF      * sqrt(M_FIVE/M_SEVEN)            * params%k(1, 2, 2)
    params%k(1, 1, 3) =  M_ONE/M_SIX * sqrt(CNST(35.0)/CNST(11.0))     * params%k(1, 3, 3)
    params%k(1, 2, 3) = -M_ONE/M_SIX * (CNST(14.0) / sqrt(CNST(11.0))) * params%k(1, 3, 3)
    params%k(2, 1, 2) = -M_HALF      * sqrt(M_SEVEN/M_NINE)            * params%k(2, 2, 2)
    params%k(2, 1, 3) =  M_HALF      * sqrt(CNST(63.0)/CNST(143.0))    * params%k(2, 3, 3)
    params%k(2, 2, 3) = -M_HALF      * (CNST(18.0)/sqrt(CNST(143.0)))  * params%k(2, 3, 3)


    ! Parameters are symmetric.
    do k = 0, 3
       do i = 1, 3
          do j = i + 1, 3
             params%h(k, j, i) = params%h(k, i, j)
             params%k(k, j, i) = params%k(k, i, j)
          enddo
       enddo
    enddo

    load_params = 0
    call pop_sub()
  end function load_params

  subroutine get_cutoff_radii(psp)
    type(hgh_type), intent(inout)     :: psp

    integer  :: ir, l, i
    FLOAT :: dincv, tmp
    FLOAT, parameter :: threshold = CNST(1.0e-4)

    call push_sub('hgh.get_cutoff_radii_psp')

    do l = 0, psp%l_max
       tmp = M_ZERO
       do i = 1, 3
          do ir = psp%g%nrval, 2, -1
             dincv = abs(psp%kb(ir, l, i))
             if(dincv > threshold) exit
          enddo
          tmp = psp%g%rofi(ir + 1)
          psp%kbr(l) = max(tmp, psp%kbr(l))
       enddo
    end do

    call pop_sub()
  end subroutine get_cutoff_radii

  ! Local pseudopotential, both in real and reciprocal space.
  function vlocalr_scalar(r, p)
    type(hgh_type), intent(in)     :: p
    FLOAT, intent(in)           :: r
    FLOAT                       :: vlocalr_scalar

    FLOAT :: r1, r2, r4, r6

    r1 = r/p%rlocal; r2 = r1**2; r4 = r2**2; r6 = r4*r2

    if(r < CNST(1.0e-7)) then
       vlocalr_scalar = - (M_TWO * p%z_val)/(sqrt(M_TWO*M_Pi)*p%rlocal) + p%c(1)
       return
    endif

    vlocalr_scalar = - (p%z_val/r)*loct_erf(r1/sqrt(M_TWO))   &
         + exp( -M_HALF*r2 ) *    &
         ( p%c(1) + p%c(2)*r2 + p%c(3)*r4 + p%c(4)*r6 )

  end function vlocalr_scalar

  function vlocalr_vector(r, p)
    type(hgh_type), intent(in)      :: p
    FLOAT, intent(in)            :: r(:)
    FLOAT, pointer               :: vlocalr_vector(:)

    integer :: i

    allocate(vlocalr_vector(size(r)))
    do i = 1, size(r)
       vlocalr_vector(i) = vlocalr_scalar(r(i), p)
    enddo

  end function vlocalr_vector

  function vlocalg(g, p)
    type(hgh_type), intent(in)     :: p
    FLOAT, intent(in)           :: g
    FLOAT                       :: vlocalg

    FLOAT :: g1, g2, g4, g6

    g1 = g*p%rlocal; g2 = g1*g1; g4 = g2*g2; g6 = g4*g2

    vlocalg = -(M_FOUR*M_Pi*p%z_val/g**2) * exp( -g2/M_TWO) +            &
         sqrt(M_EIGHT*M_Pi**3) * p%rlocal**3 * exp( -g2/M_TWO) *    &
         ( p%c(1) + p%c(2)*(M_THREE - g2) + p%c(3)*(CNST(15.0) - M_TEN*g2 + g4) + &
         p%c(4)*(CNST(105.0) -CNST(105.0)*g2 + CNST(21.0)*g4 - g6) )

  end function vlocalg

  function projectorr_scalar(r, p, i, l)
    type(hgh_type), intent(in)     :: p
    FLOAT, intent(in)           :: r
    integer, intent(in)            :: i, l
    FLOAT                       :: projectorr_scalar

    FLOAT :: x, y, rr

    x = l + real(4*i-1, PRECISION)/M_TWO
    y = loct_gamma(x); x = sqrt(y)
    if(l==0 .and. i==1) then
       rr = M_ONE
    else
       rr = r ** (l + 2*(i-1))
    endif

    projectorr_scalar = sqrt(M_TWO) * rr * exp(-r**2/(M_TWO*p%rc(l)**2)) / &
         (  p%rc(l)**(l + real(4*i-1, PRECISION)/M_TWO) * x )

  end function projectorr_scalar

  function projectorr_vector(r, p, i, l)
    type(hgh_type), intent(in)     :: p
    FLOAT, intent(in)           :: r(:)
    integer, intent(in)            :: i, l
    FLOAT, pointer              :: projectorr_vector(:)

    integer :: j

    allocate(projectorr_vector(size(r)))
    do j=1, size(r)
       projectorr_vector(j) = projectorr_scalar(r(j), p, i, l)
    enddo

  end function projectorr_vector

  function projectorg(g, p, i, l)
    type(hgh_type), intent(in)     :: p
    FLOAT, intent(in)           :: g
    integer, intent(in)            :: i, l
    FLOAT                       :: projectorg

    !FLOAT, external :: gamma
    FLOAT :: pif, ex

    pif = M_Pi**(M_FIVE/M_FOUR)

    ex = exp( M_HALF*(g*p%rc(l))**2 )

    select case(l)
    case(0)
       select case(i)
       case(1)
          projectorg = ( M_FOUR*sqrt(M_TWO*p%rc(0)**3)*pif ) / ex
       case(2)
          projectorg = ( sqrt(M_EIGHT*2*p%rc(0)**3/CNST(15.0))*pif * &
               (M_THREE - (g*p%rc(0))**2) ) / ex
       case(3)
          projectorg = ( CNST(16.0)*sqrt(M_TWO*p%rc(0)**3/CNST(105.0)) * pif * &
               (CNST(15.0) - M_TEN*g**2*p%rc(0)**2 + g**4*p%rc(0)**2) ) / (M_THREE*ex)
       end select

    case(1)
       select case(i)
       case(1)
          projectorg = ( M_EIGHT*sqrt(p%rc(1)**5/M_THREE)*pif*g ) / ex
       case(2)
          projectorg = ( CNST(16.0)*sqrt(p%rc(1)**5/CNST(105.0))* pif * g * &
               ( M_FIVE - (g*p%rc(1))**2 ) ) / ex
       case(3)
          projectorg = ( CNST(32.0)*sqrt(p%rc(1)**5/CNST(1155.0))* pif * g * &
               ( CNST(35.0) - CNST(14.0)*g**2*p%rc(1)**2 + (g*p%rc(1))**4 ) ) / &
               (M_THREE*ex)
       end select

    case(2)
       select case(i)
       case(1)
          projectorg = ( M_EIGHT * sqrt(M_TWO*p%rc(2)**7/CNST(15.0)) * pif * g**2 ) / ex
       case(2)
          projectorg = ( CNST(16.0) * sqrt(M_TWO*p%rc(2)**7/CNST(105.0)) * pif * g**2 * &
               (M_SEVEN - g**2*p%rc(2)**2) ) / (M_THREE*ex)
       case(3)
          projectorg = M_ZERO ! ??
       end select

    case(3)
       ! This should be checked. Probably will not be needed in an near future...
    end select

  end function projectorg

  subroutine solve_schroedinger(psp)
    type(hgh_type), intent(inout)     :: psp

    integer :: iter, ir, l, nnode, nprin, i, j, irr, n, k
    FLOAT :: vtot, a2b4, diff, nonl
    FLOAT, allocatable :: prev(:, :), rho(:, :), ve(:, :)
    FLOAT, parameter :: tol = CNST(1.0e-4)
    DOUBLE :: e, z, dr, rmax
    DOUBLE, allocatable :: s(:), hato(:), g(:), y(:)

    call push_sub('hgh.solve_schroedinger')

    ! Allocations...
    allocate(s(psp%g%nrval), hato(psp%g%nrval), g(psp%g%nrval), y(psp%g%nrval), prev(psp%g%nrval, 1), &
         rho(psp%g%nrval, 1), ve(psp%g%nrval, 1))
    hato = M_ZERO; g = M_ZERO;  y = M_ZERO; rho = M_ZERO; ve = M_ZERO

    ! These numerical parameters have to be fixed for egofv to work.
    s(2:psp%g%nrval) = psp%g%drdi(2:psp%g%nrval)*psp%g%drdi(2:psp%g%nrval)
    s(1) = s(2)
    a2b4 = M_FOURTH*psp%g%a**2

    ! Let us be a bit informative.
    if(conf%verbose > 20 .and. mpiv%node == 0) then
       message(1) = '      Calculating atomic pseudo-eigenfunctions for specie ' // psp%atom_name // '....'
       call write_info(1)
    endif

    ! "Double" self consistent loop: nonlocal and xc parts have to be calculated
    ! self-consistently.
    diff = CNST(1e5)
    iter = 0
    self_consistent: do while( diff > tol )
       prev = rho
       iter = iter + 1
       do n = 1, psp%conf%p
          l = psp%conf%l(n)
          do ir = 2, psp%g%nrval
             vtot = 2*psp%vlocal(ir) + ve(ir, 1) + dble(l*(l + 1))/(psp%g%rofi(ir)**2)
             nonl = M_ZERO
             if(iter>1 .and. psp%l_max >=0 .and. psp%rphi(ir, n) > 1.0e-7) then
                do i = 1, 3
                   do j = 1, 3
                      do irr = 2, psp%g%nrval
                         nonl = nonl + psp%h(l, i, j)*psp%kb(ir, l, i)* &
                              psp%g%drdi(irr)*psp%g%rofi(irr)*psp%rphi(irr, n)*psp%kb(irr,l,j)
                      enddo
                   enddo
                enddo
                nonl = 2*nonl/psp%rphi(ir, n)*psp%g%rofi(ir)
             endif
             vtot = vtot + nonl
             hato(ir) = vtot*s(ir) + a2b4
          end do
          hato(1) = hato(2)
          ! We will assume there is only the possibility of two equal l numbers...
          nnode = 1
          do k = 1, n - 1
             if(psp%conf%l(k)==psp%conf%l(n)) nnode = 2
          enddo
          nprin = l + 1
          if(iter == 1) then
             e = -((psp%z_val/dble(nprin))**2); z = psp%z_val
          else
             e = psp%eigen(n); z = psp%z_val
          endif
          dr = -CNST(1.0e5); rmax = psp%g%rofi(psp%g%nrval)
          call egofv(hato, s, psp%g%nrval, e, g, y, l, z, &
               real(psp%g%a, 8), real(psp%g%b, 8), rmax, nprin, nnode, dr)
          psp%eigen(n) = e

          psp%rphi(2:psp%g%nrval, n) = g(2:psp%g%nrval) * sqrt(psp%g%drdi(2:psp%g%nrval))
          psp%rphi(1, n) = psp%rphi(2, n)
       end do
       rho = M_ZERO
       do n = 1, psp%conf%p
          rho(1:psp%g%nrval, 1) = rho(1:psp%g%nrval, 1) + psp%conf%occ(n,1)*psp%rphi(1:psp%g%nrval, n)**2
       enddo
       if(iter>1) rho = M_HALF*(rho + prev)
       diff = sqrt(sum(psp%g%drdi(2:psp%g%nrval)*(rho(2:psp%g%nrval, 1)-prev(2:psp%g%nrval, 1))**2))
       !if(conf%verbose>20 .and. mpiv%node == 0) then
       !  write(message(1),'(a,i4,a,e10.2)') '      Iter =',iter,'; Diff =',diff
       !  call write_info(1)
       !endifq
       call atomhxc('LDA', psp%g, 1, rho, ve)

    enddo self_consistent

    if(conf%verbose > 20 .and. mpiv%node == 0) then
       message(1) = '      Done.'
       call write_info(1)
    endif

    !  checking normalization of the calculated wave functions
    !do l = 0, psp%l_max_occ
    do n = 1, psp%conf%p
       e = sqrt(sum(psp%g%drdi(2:psp%g%nrval)*psp%rphi(2:psp%g%nrval, n)**2))
       e = abs(e - M_ONE)
       if (e > CNST(1.0e-5) .and. conf%verbose > 0) then
          write(message(1), '(a,i2,a)') "Eigenstate for n = ", n , ' is not normalized'
          write(message(2), '(a, f12.6,a)') '(abs(1-norm) = ', e, ')'
          call write_warning(2)
       end if
    end do

    ! Output in Ha and not in stupid Rydbergs.
    psp%eigen = psp%eigen / M_TWO

    ! Deallocations.
    deallocate(s, hato, g, y, rho, prev)

    call pop_sub()
  end subroutine solve_schroedinger

  subroutine hgh_debug(psp, dir)
    type(hgh_type), intent(in) :: psp
    character(len=*), intent(in) :: dir

    integer :: hgh_unit, loc_unit, dat_unit, kbp_unit, wav_unit, i, l, k
    character(len=256) :: dirname

    call push_sub('hgh.hgh_debug')

    ! Open files.
    dirname = trim(dir)//'/hgh.'//trim(psp%atom_name)
    call io_mkdir(trim(dir))
    hgh_unit = io_open(trim(dirname)//'/hgh', action='write')
    loc_unit = io_open(trim(dirname)//'/local', action='write')
    dat_unit = io_open(trim(dirname)//'/info', action='write')
    kbp_unit = io_open(trim(dirname)//'/nonlocal', action='write')
    wav_unit = io_open(trim(dirname)//'/wave', action='write')

    ! Writes down the input file, to be checked agains SHARE_OCTOPUS/PP/HGH/ATOM_NAME.hgh
    write(hgh_unit,'(a5,i6,5f12.6)') psp%atom_name, psp%z_val, psp%rlocal, psp%c(1:4)
    write(hgh_unit,'(  11x,4f12.6)') psp%rc(0), (psp%h(0,i,i), i = 1, 3)
    do k = 1, 3
       write(hgh_unit,'(  11x,4f12.6)') psp%rc(k), (psp%h(k, i, i), i = 1, 3)
       write(hgh_unit,'(  23x,4f12.6)')            (psp%k(k, i, i), i = 1, 3)
    enddo

    ! Writes down some info.
    write(dat_unit,'(a,i3)')        'lmax  = ', psp%l_max
    if(psp%l_max >= 0) then
       write(dat_unit,'(a,4f14.6)') 'kbr   = ', psp%kbr
    endif
    write(dat_unit,'(a,5f14.6)')    'eigen = ', psp%eigen
    write(dat_unit,'(a,5f14.6)')    'occ   = ', psp%conf%occ(1:psp%conf%p, 1)
    ! Writes down local part.
    do i = 1, psp%g%nrval
       write(loc_unit, *) psp%g%rofi(i), psp%vlocal(i)
    enddo

    ! Writes down nonlocal part.
    if(psp%l_max >=0) then
       do i = 1, psp%g%nrval
          write(kbp_unit, '(10es14.4)') psp%g%rofi(i), ( (psp%kb(i, l, k) ,k = 1, 3),l = 0, psp%l_max)
       enddo
    endif

    ! And the pseudo-wavefunctions.
    do i = 1, psp%g%nrval
       write(wav_unit, *) psp%g%rofi(i), (psp%rphi(i, l), l = 1, psp%conf%p)
    enddo

    ! Closes files and exits
    call io_close(hgh_unit); call io_close(loc_unit); call io_close(wav_unit)
    call io_close(dat_unit); call io_close(kbp_unit)

    call pop_sub()
  end subroutine hgh_debug

end module hgh
