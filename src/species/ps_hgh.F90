!! Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
!! Copyright (C) 2020 N. Tancogne-Dejean
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

#include "global.h"

module ps_hgh_oct_m
  !< For information about the Hartwinger-Goedecker-Hutter pseudopotentials, take a look at:
  !!  (1) S. Goedecker, M. Teter and J. Hutter, Phys. Rev. B 54, 1703 (1996).
  !!  (2) C. Hartwinger, S. Goedecker and J. Hutter, Phys. Rev. B 58, 3641 (1998).
  use atomic_oct_m
  use global_oct_m
  use io_oct_m
  use loct_math_oct_m
  use logrid_oct_m
  use messages_oct_m
  use namespace_oct_m
  use profiling_oct_m

  implicit none

  private
  public  ::      &
    ps_hgh_t,     &
    hgh_init,     &
    hgh_process,  &
    hgh_debug,    &
    hgh_end,      &
    hgh_get_eigen

  !> The following data type contains:
  !!   (a) the pseudopotential parameters, as read from a *.hgh file,
  !!   (b) auxiliary intermediate functions, to store stuff before passing it to the "ps" variable.
  type ps_hgh_t
    ! Components are public by default

    character(len=256), private :: title
    integer,            private :: pspdat ! date of creation of PP (DDMMYY)
    integer,            private :: pspcod ! code for the pseudopotential (3 for .hgh)
    integer,            private :: pspxc  ! exchange-correlation used to generate the psp
    integer                     :: lmax   ! Maximum l to use
    integer,            private :: mmax   ! Maximum number of points in real space grid
    FLOAT,              private :: r2well ! ??

    !< HGH parameters.
    character(len=5), private :: atom_name
    integer          :: z_val
    FLOAT, private   :: rlocal
    FLOAT, private   :: rc(0:3)
    FLOAT, private   :: c(1:4)
    FLOAT            :: h(0:3, 1:3, 1:3)
    FLOAT            :: k(0:3, 1:3, 1:3)

    type(valconf_t)  :: conf
    integer          :: l_max     !< Maximum l for the Kleinman-Bylander component.

    FLOAT, allocatable   :: vlocal(:) !< Local potential
    FLOAT, allocatable   :: kb(:,:,:) !< KB projectors
    FLOAT, allocatable   :: kbr(:)    !< KB radii
    FLOAT, allocatable   :: rphi(:,:)
    FLOAT, allocatable, private :: eigen(:)

    !> Logarithmic grid parameters
    type(logrid_t) :: g
  end type ps_hgh_t

  FLOAT, parameter :: eps = CNST(1.0e-8)

contains

  ! ---------------------------------------------------------
  subroutine hgh_init(psp, filename, namespace)
    type(ps_hgh_t),    intent(inout) :: psp
    character(len=*),  intent(in)    :: filename
    type(namespace_t), intent(in)    :: namespace

    integer :: iunit, i

    PUSH_SUB(hgh_init)

    iunit = io_open(trim(filename), action='read', form='formatted', status='old')
    i = load_params(iunit, psp, namespace)
    if(i /= 0) then
      call messages_write('Error reading hgh file')
      call messages_fatal(namespace=namespace)
    end if
    call io_close(iunit)

    ! Finds out psp%l_max. The most special cases are H, He, Li_sc and Be_sc, where psp%l_max = -1.
    psp%l_max = 0
    do while(psp%rc(psp%l_max) > CNST(0.01))
      psp%l_max = psp%l_max + 1
      if(psp%l_max > 3) exit
    end do
    psp%l_max = psp%l_max - 1

    ! Initializes the logarithmic grid. Parameters are hard-coded.
    call logrid_init(psp%g, LOGRID_PSF, CNST(3.0e-2), CNST(4.0e-4), 431)

    ! Allocation of stuff.
    SAFE_ALLOCATE(psp%vlocal(1:psp%g%nrval))
    psp%vlocal = M_ZERO
    if(psp%l_max >= 0) then
      SAFE_ALLOCATE(psp%kbr(0:psp%l_max))
      SAFE_ALLOCATE(psp%kb(1:psp%g%nrval, 0:psp%l_max, 1:3))
      psp%kbr = M_ZERO
      psp%kb = M_ZERO
    end if

    POP_SUB(hgh_init)
  end subroutine hgh_init


  ! ---------------------------------------------------------
  subroutine hgh_end(psp)
    type(ps_hgh_t), intent(inout) :: psp

    PUSH_SUB(hgh_end)

    if(psp%l_max >= 0) then
      SAFE_DEALLOCATE_A(psp%kbr)
      SAFE_DEALLOCATE_A(psp%kb)
    end if
    SAFE_DEALLOCATE_A(psp%vlocal)
    SAFE_DEALLOCATE_A(psp%rphi)
    SAFE_DEALLOCATE_A(psp%eigen)
    call logrid_end(psp%g)

    POP_SUB(hgh_end)
  end subroutine hgh_end


  ! ---------------------------------------------------------
  subroutine hgh_process(psp, namespace)
    type(ps_hgh_t),    intent(inout) :: psp
    type(namespace_t), intent(in)    :: namespace

    integer :: l, i, ierr

    PUSH_SUB(hgh_process)

    SAFE_ALLOCATE(psp%rphi(1:psp%g%nrval, 1:psp%conf%p))
    SAFE_ALLOCATE(psp%eigen(1:psp%conf%p))
    psp%rphi = M_ZERO
    psp%eigen = M_ZERO


    ! Fixes the local potential
    call vlocalr_scalar(psp%g%rofi, psp%g%nrval, psp, psp%vlocal)

    ! And the projectors
    do l = 0, psp%l_max
      do i = 1, 3
        call projectorr_scalar(psp%g%rofi, psp%g%nrval, psp, i, l, psp%kb(:, l, i))
      end do
    end do

    ! get the pseudoatomic eigenfunctions (WARNING: This is not correctly done yet: "some" wavefunctions
    ! are obtained, but not the real ones!!!
    call solve_schroedinger(psp, ierr, namespace)
    if(ierr /= 0) then ! If the wavefunctions could not be found, we set its number to zero.
      write(message(1),'(a)') 'The algorithm that calculates atomic wavefunctions could not'
      write(message(2),'(a)') 'do its job. The program will continue, but expect poor'
      write(message(3),'(a)') 'convergence properties.'
      call messages_warning(3, namespace=namespace)
      psp%conf%p = 0
    end if

    ! Define the KB-projector cut-off radii
    call get_cutoff_radii(psp)

    POP_SUB(hgh_process)
  end subroutine hgh_process

  ! ---------------------------------------------------------
  subroutine hgh_get_eigen(psp, eigen)
    type(ps_hgh_t), intent(in)  :: psp
    FLOAT,          intent(out) :: eigen(:,:)

    integer :: i
    
    PUSH_SUB(hgh_get_eigen)

    do i = 1, psp%conf%p
      eigen(i, :) = psp%eigen(i)
    end do

    POP_SUB(hgh_get_eigen)
  end subroutine hgh_get_eigen
  
  ! ---------------------------------------------------------
  function load_params(unit, params, namespace)
    integer,             intent(in)  :: unit        ! where to read from
    type(ps_hgh_t),      intent(out) :: params      ! obvious
    type(namespace_t),   intent(in)  :: namespace

    integer                     :: load_params ! 0 if success,
    ! 1 otherwise.

    integer :: i, iostat, j, k, nn, nnonloc, lloc
    character(len=VALCONF_STRING_LENGTH) :: line

    PUSH_SUB(load_params)

    ! Set initially everything to zero.
    params%c(1:4) = M_ZERO
    params%rlocal = M_ZERO
    params%rc = M_ZERO
    params%h = M_ZERO
    params%k = M_ZERO

    read(unit, *) params%title
    read(unit, *) params%conf%z, params%z_val, params%pspdat
    read(unit, *, iostat=iostat) params%pspcod, params%pspxc, params%l_max, lloc, params%mmax, params%r2well    

    select case(params%pspcod)
    case(3) 

      ! Reads the file in a hopefully smart way
      iostat = 1
      j = 5
      read(unit,'(a)') line
      do while((iostat /= 0) .and. (j > 0))
        j = j - 1
        read(line, *, iostat=iostat) params%rlocal, params%c(1:j)
      end do
      if(j<1) read(line, *, iostat=iostat) params%atom_name, params%z_val, params%rlocal
      if( iostat /= 0 ) then
        load_params = 1
        POP_SUB(load_params)
        return
      end if

      read(unit,'(a)', iostat = iostat) line
      if(iostat /= 0) then
        load_params = 0
        POP_SUB(load_params)
        return
      end if
      iostat = 1
      j = 4
      do while((iostat /= 0) .and. (j > 0))
        j = j - 1
        read(line, *, iostat=iostat) params%rc(0), (params%h(0, i, i), i = 1, j)
      end do
      if(j < 0) then
        load_params = 2
        POP_SUB(load_params)
        return
      end if

      kloop: do k = 1, 3
        read(unit, '(a)', iostat = iostat) line
        if(iostat /= 0) exit kloop
        iostat = 1
        j = 4
        do while((iostat /= 0) .and. (j > 0))
          j = j - 1
          read(line, *, iostat = iostat) params%rc(k), (params%h(k, i, i), i = 1, j)
        end do
        if(abs(params%rc(k)) <= M_EPSILON) exit kloop
        read(unit, '(a)') line
        iostat = 1
        j = 4
        do while((iostat /= 0) .and. (j>0))
          j = j - 1
          read(line, *, iostat = iostat) (params%k(k, i, i), i = 1, 3)
        end do
      end do kloop

      ! Fill in the rest of the parameter matrices...
      params%h(0, 1, 2) = -M_HALF      * sqrt(M_THREE/M_FIVE)            * params%h(0, 2, 2)
      params%h(0, 1, 3) =  M_HALF      * sqrt(M_FIVE/CNST(21.0))         * params%h(0, 3, 3)
      params%h(0, 2, 3) = -M_HALF      * sqrt(CNST(100.0)/CNST(63.0))    * params%h(0, 3, 3)
      params%h(1, 1, 2) = -M_HALF      * sqrt(M_FIVE/CNST(7.0))            * params%h(1, 2, 2)
      params%h(1, 1, 3) =  M_ONE/CNST(6.0) * sqrt(CNST(35.0)/CNST(11.0))     * params%h(1, 3, 3)
      params%h(1, 2, 3) = -M_ONE/CNST(6.0) * (CNST(14.0) / sqrt(CNST(11.0))) * params%h(1, 3, 3)
      params%h(2, 1, 2) = -M_HALF      * sqrt(CNST(7.0)/CNST(9.0))            * params%h(2, 2, 2)
      params%h(2, 1, 3) =  M_HALF      * sqrt(CNST(63.0)/CNST(143.0))    * params%h(2, 3, 3)
      params%h(2, 2, 3) = -M_HALF      * (CNST(18.0)/sqrt(CNST(143.0)))  * params%h(2, 3, 3)

      params%k(0, 1, 2) = -M_HALF      * sqrt(M_THREE/M_FIVE)            * params%k(0, 2, 2)
      params%k(0, 1, 3) =  M_HALF      * sqrt(M_FIVE/CNST(21.0))         * params%k(0, 3, 3)
      params%k(0, 2, 3) = -M_HALF      * sqrt(CNST(100.0)/CNST(63.0))    * params%k(0, 3, 3)
      params%k(1, 1, 2) = -M_HALF      * sqrt(M_FIVE/CNST(7.0))            * params%k(1, 2, 2)
      params%k(1, 1, 3) =  M_ONE/CNST(6.0) * sqrt(CNST(35.0)/CNST(11.0))     * params%k(1, 3, 3)
      params%k(1, 2, 3) = -M_ONE/CNST(6.0) * (CNST(14.0) / sqrt(CNST(11.0))) * params%k(1, 3, 3)
      params%k(2, 1, 2) = -M_HALF      * sqrt(CNST(7.0)/CNST(9.0))            * params%k(2, 2, 2)
      params%k(2, 1, 3) =  M_HALF      * sqrt(CNST(63.0)/CNST(143.0))    * params%k(2, 3, 3)
      params%k(2, 2, 3) = -M_HALF      * (CNST(18.0)/sqrt(CNST(143.0)))  * params%k(2, 3, 3)


    case(10)

      read(unit, *) params%rlocal, nn, params%c(1)
      read(unit, *) nnonloc
      read(unit, '(a)', iostat = iostat) line
      kloop2 : do k = 0, 3
        if(iostat /= 0) exit kloop2
        iostat = 1
        j = 4
        do while((iostat /= 0) .and. (j > 0))
          j = j - 1
          read(line, *, iostat = iostat) params%rc(k), nn, (params%h(k, 1, i), i = 1, j)
        end do
        if(abs(params%rc(k)) <= M_EPSILON) exit kloop2
        read(unit, '(a)') line
        select case(nn)
        case(3)
          read(line, *) (params%h(1, 2, i), i = 2, 3)
          read(unit, '(a)') line
          read(line, *) params%h(1, 3, 3)
          !Reading the k matrix, if possible
          read(unit, '(a)', iostat=iostat) line
          if(iostat /= 0) exit kloop2
          read(line, *, iostat=iostat) (params%k(1, 1, i), i = 1, 3)
          if(iostat /= 0) continue !No k matrix 
          read(unit, '(a)') line
          read(line, *) (params%k(1, 2, i), i = 2, 3)
          read(unit, '(a)') line
          read(line, *) params%k(1, 3, 3)
          read(unit, '(a)', iostat = iostat) line
        case(2)
          read(line, *) params%h(0, 2, 2)
          !Reading the k matrix, if possible
          read(unit, '(a)', iostat=iostat) line
          if(iostat /= 0) exit kloop2
          read(line, *, iostat=iostat) (params%k(1, 1, i), i = 1, 2)
          if(iostat /= 0) continue !No k matrix 
          read(unit, '(a)') line
          read(line, *) params%k(1, 2, 2)   
          read(unit, '(a)', iostat = iostat) line
        case default
          message(1) = "Error parsing the pseudopotential"
          call messages_fatal(1, namespace=namespace)
        end select
      end do kloop2

    case default
      message(1) = "Inconsistency in pseudopotential file:"
      write(message(2),'(a,i2)') "  expecting pspcod = 3, but found ", params%pspcod
      call messages_fatal(2, namespace=namespace)
    end select


    ! Parameters are symmetric.
    do k = 0, 3
      do i = 1, 3
        do j = i + 1, 3
          params%h(k, j, i) = params%h(k, i, j)
          params%k(k, j, i) = params%k(k, i, j)
        end do
      end do
    end do

    load_params = 0
    POP_SUB(load_params)
  end function load_params


  ! ---------------------------------------------------------
  subroutine get_cutoff_radii(psp)
    type(ps_hgh_t), intent(inout)     :: psp

    integer  :: ir, l, i
    FLOAT :: dincv, tmp
    FLOAT, parameter :: threshold = CNST(1.0e-4)

    PUSH_SUB(get_cutoff_radii)

    do l = 0, psp%l_max
      tmp = M_ZERO
      do i = 1, 3
        do ir = psp%g%nrval, 2, -1
          dincv = abs(psp%kb(ir, l, i))
          if(dincv > threshold) exit
        end do
        tmp = psp%g%rofi(ir + 1)
        psp%kbr(l) = max(tmp, psp%kbr(l))
      end do
    end do

    POP_SUB(get_cutoff_radii)
  end subroutine get_cutoff_radii


  ! ---------------------------------------------------------
  ! Local pseudopotential, both in real and reciprocal space.
  ! See Eq. (1)
  subroutine vlocalr_scalar(r, np, p, vloc)
    type(ps_hgh_t), intent(in)    :: p
    FLOAT,          intent(in)    :: r(:)
    integer,        intent(in)    :: np
    FLOAT,          intent(inout) :: vloc(:)

    integer :: ip
    FLOAT :: r1, r2, r4, r6

    PUSH_SUB(vlocalr_scalar)

    do ip = 1, np
      if(r(ip) < CNST(1.0e-7)) then
        vloc(ip) = - (M_TWO * p%z_val)/(sqrt(M_TWO*M_Pi)*p%rlocal) + p%c(1)
      else
        r1 = r(ip)/p%rlocal
        r2 = r1**2
        r4 = r2**2
        r6 = r4*r2

        vloc(ip) = - (p%z_val/r(ip))*loct_erf(r1/sqrt(M_TWO))   &
          + exp( -M_HALF*r2 ) *    &
          ( p%c(1) + p%c(2)*r2 + p%c(3)*r4 + p%c(4)*r6 )

      end if
    end do

    POP_SUB(vlocalr_scalar)
  end subroutine vlocalr_scalar


  ! ---------------------------------------------------------
  ! Projector in real space, see Eq. 3
  subroutine projectorr_scalar(r, np, p, i, l, proj)
    type(ps_hgh_t), intent(in)    :: p
    FLOAT,          intent(in)    :: r(:)
    integer,        intent(in)    :: np
    integer,        intent(in)    :: i
    integer,        intent(in)    :: l
    FLOAT,          intent(inout) :: proj(:)

    integer :: ip
    FLOAT :: x, y, rr

    PUSH_SUB(projectorr_scalar)

    x = l + TOFLOAT(4*i-1)/M_TWO
    y = loct_gamma(x)
    x = sqrt(y)
    rr = M_ONE

    do ip = 1, np

      if(l /=0 .or. i /= 1) then
        rr = r(ip) ** (l + 2*(i-1))
      end if

      proj(ip) = sqrt(M_TWO) * rr * exp(-r(ip)**2/(M_TWO*p%rc(l)**2)) / &
            (  p%rc(l)**(l + TOFLOAT(4*i-1)/M_TWO) * x )

    end do

    POP_SUB(projectorr_scalar)
  end subroutine projectorr_scalar


  ! ---------------------------------------------------------
  subroutine solve_schroedinger(psp, ierr, namespace)
    type(ps_hgh_t),    intent(inout) :: psp
    integer,           intent(out)   :: ierr
    type(namespace_t), intent(in)    :: namespace

    integer :: iter, ir, l, nnode, nprin, i, j, irr, n, k
    FLOAT :: vtot, a2b4, diff, nonl
    FLOAT, allocatable :: prev(:, :), rho(:, :), ve(:, :)
    FLOAT, parameter :: tol = CNST(1.0e-4)
    REAL_DOUBLE :: e, z, dr, rmax
    REAL_DOUBLE, allocatable :: s(:), hato(:), g(:), y(:)

    PUSH_SUB(solve_schroedinger)

    ierr = 0

    ! Allocations...
    SAFE_ALLOCATE(   s(1:psp%g%nrval))
    SAFE_ALLOCATE(hato(1:psp%g%nrval))
    SAFE_ALLOCATE(   g(1:psp%g%nrval))
    SAFE_ALLOCATE(   y(1:psp%g%nrval))
    SAFE_ALLOCATE(prev(1:psp%g%nrval, 1:1))
    SAFE_ALLOCATE( rho(1:psp%g%nrval, 1:1))
    SAFE_ALLOCATE(  ve(1:psp%g%nrval, 1:1))
    hato = M_ZERO
    g = M_ZERO
    y = M_ZERO
    rho = M_ZERO
    ve = M_ZERO

    ! These numerical parameters have to be fixed for egofv to work.
    s(2:psp%g%nrval) = psp%g%drdi(2:psp%g%nrval)*psp%g%drdi(2:psp%g%nrval)
    s(1) = s(2)
    a2b4 = M_FOURTH*psp%g%a**2
    
    ! "Double" self-consistent loop: nonlocal and XC parts have to be calculated
    ! self-consistently.
    diff = CNST(1e5)
    iter = 0
    self_consistent: do while( diff > tol )
      prev = rho
      iter = iter + 1
      do n = 1, psp%conf%p
        l = psp%conf%l(n)
        do ir = 2, psp%g%nrval
          vtot = 2*psp%vlocal(ir) + ve(ir, 1) + TOFLOAT(l*(l + 1))/(psp%g%rofi(ir)**2)
          nonl = M_ZERO
          if(iter>2 .and. psp%l_max >=0 .and. psp%rphi(ir, n) > CNST(1.0e-7)) then
            do i = 1, 3
              do j = 1, 3
                do irr = 2, psp%g%nrval
                  nonl = nonl + psp%h(l, i, j)*psp%kb(ir, l, i)* &
                    psp%g%drdi(irr)*psp%g%rofi(irr)*psp%rphi(irr, n)*psp%kb(irr,l,j)
                end do
              end do
            end do
            nonl = 2*nonl/psp%rphi(ir, n)*psp%g%rofi(ir)
          end if
          vtot = vtot + nonl
          hato(ir) = vtot*s(ir) + a2b4
        end do
        hato(1) = hato(2)
        ! We will assume there is only the possibility of two equal l numbers...
        nnode = 1
        do k = 1, n - 1
          if(psp%conf%l(k)==psp%conf%l(n)) nnode = 2
        end do
        nprin = l + 1
        if(iter == 1) then
          e = -((psp%z_val/TOFLOAT(nprin))**2)
          z = psp%z_val
        else
          e = psp%eigen(n)
          z = psp%z_val
        end if
        dr = -CNST(1.0e5)
        rmax = psp%g%rofi(psp%g%nrval)
        call egofv(hato, s, psp%g%nrval, e, g, y, l, z, &
          real(psp%g%a, 8), real(psp%g%b, 8), rmax, nprin, nnode, dr, ierr)
        if(ierr /= 0) exit self_consistent
        psp%eigen(n) = e

        psp%rphi(2:psp%g%nrval, n) = g(2:psp%g%nrval) * sqrt(psp%g%drdi(2:psp%g%nrval))
        psp%rphi(1, n) = psp%rphi(2, n)
      end do
      rho = M_ZERO
      do n = 1, psp%conf%p
        rho(1:psp%g%nrval, 1) = rho(1:psp%g%nrval, 1) + psp%conf%occ(n,1)*psp%rphi(1:psp%g%nrval, n)**2
      end do
      if(iter>1) rho = M_HALF*(rho + prev)
      diff = sqrt(sum(psp%g%drdi(2:psp%g%nrval)*(rho(2:psp%g%nrval, 1)-prev(2:psp%g%nrval, 1))**2))
      call atomhxc('LDA', psp%g, 1, rho, ve)

    end do self_consistent

    if(ierr  ==  0) then
      !  checking normalization of the calculated wavefunctions
      !do l = 0, psp%l_max_occ
      do n = 1, psp%conf%p
        e = sqrt(sum(psp%g%drdi(2:psp%g%nrval)*psp%rphi(2:psp%g%nrval, n)**2))
        e = abs(e - M_ONE)
        if (e > CNST(1.0e-5)) then
          write(message(1), '(a,i2,a)') "Eigenstate for n = ", n , ' is not normalized'
          write(message(2), '(a, f12.6,a)') '(abs(1-norm) = ', e, ')'
          call messages_warning(2, namespace=namespace)
        end if
      end do

      ! Output in Ha and not in stupid Rydbergs.
      psp%eigen = psp%eigen / M_TWO
    end if

    ! Deallocations.
    SAFE_DEALLOCATE_A(s)
    SAFE_DEALLOCATE_A(hato)
    SAFE_DEALLOCATE_A(g)
    SAFE_DEALLOCATE_A(y)
    SAFE_DEALLOCATE_A(prev)
    SAFE_DEALLOCATE_A(rho)
    SAFE_DEALLOCATE_A(ve)

    POP_SUB(solve_schroedinger)
  end subroutine solve_schroedinger


  ! ---------------------------------------------------------
  subroutine hgh_debug(psp, dir, namespace)
    type(ps_hgh_t),    intent(in) :: psp
    character(len=*),  intent(in) :: dir
    type(namespace_t), intent(in) :: namespace

    integer :: hgh_unit, loc_unit, dat_unit, kbp_unit, wav_unit, i, l, k
    character(len=256) :: dirname

    PUSH_SUB(hgh_debug)

    ! Open files.
    dirname = trim(dir)//'/hgh.'//trim(psp%atom_name)
    call io_mkdir(trim(dirname), namespace)
    hgh_unit = io_open(trim(dirname)//'/hgh', namespace, action='write')
    loc_unit = io_open(trim(dirname)//'/local', namespace, action='write')
    dat_unit = io_open(trim(dirname)//'/info', namespace, action='write')
    kbp_unit = io_open(trim(dirname)//'/nonlocal', namespace, action='write')
    wav_unit = io_open(trim(dirname)//'/wave', namespace, action='write')

    ! Writes down the input file, to be checked against SHARE_DIR/pseudopotentials/HGH/ATOM_NAME.hgh
    write(hgh_unit,'(a5,i6,5f12.6)') psp%atom_name, psp%z_val, psp%rlocal, psp%c(1:4)
    write(hgh_unit,'(  11x,4f12.6)') psp%rc(0), (psp%h(0,i,i), i = 1, 3)
    do k = 1, 3
      write(hgh_unit,'(  11x,4f12.6)') psp%rc(k), (psp%h(k, i, i), i = 1, 3)
      write(hgh_unit,'(  23x,4f12.6)')            (psp%k(k, i, i), i = 1, 3)
    end do

    ! Writes down some info.
    write(dat_unit,'(a,i3)')        'lmax  = ', psp%l_max
    if(psp%l_max >= 0) then
      write(dat_unit,'(a,4f14.6)') 'kbr   = ', psp%kbr
    end if
    write(dat_unit,'(a,5f14.6)')    'eigen = ', psp%eigen
    write(dat_unit,'(a,5f14.6)')    'occ   = ', psp%conf%occ(1:psp%conf%p, 1)
    ! Writes down local part.
    do i = 1, psp%g%nrval
      write(loc_unit, *) psp%g%rofi(i), psp%vlocal(i)
    end do

    ! Writes down nonlocal part.
    if(psp%l_max >=0) then
      do i = 1, psp%g%nrval
        write(kbp_unit, '(10es14.4)') psp%g%rofi(i), ( (psp%kb(i, l, k) ,k = 1, 3),l = 0, psp%l_max)
      end do
    end if

    ! And the pseudo-wavefunctions.
    do i = 1, psp%g%nrval
      write(wav_unit, *) psp%g%rofi(i), (psp%rphi(i, l), l = 1, psp%conf%p)
    end do

    ! Closes files and exits
    call io_close(hgh_unit)
    call io_close(loc_unit)
    call io_close(wav_unit)
    call io_close(dat_unit)
    call io_close(kbp_unit)

    POP_SUB(hgh_debug)
  end subroutine hgh_debug

end module ps_hgh_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
