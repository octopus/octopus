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
!! -*- coding: utf-8 mode: f90 -*-
!! $Id$

#include "global.h"

module ps_tm_m
  use global_m
  use messages_m
  use atomic_m
  use io_m
  use logrid_m
  use ps_in_grid_m
  use ps_tm_file_m

  implicit none

  private
  public ::     &
    ps_tm_t,       &
    ps_tm_init,    &
    ps_tm_end,     &
    ps_tm_process

  type ps_tm_t

    type(ps_tm_file_t) :: tm_file
    type(ps_in_grid_t) :: ps_grid

    type(valconf_t)    :: conf
    FLOAT, pointer     :: eigen(:, :)
    integer            :: ispin
  end type ps_tm_t

contains

  ! ---------------------------------------------------------
  subroutine ps_tm_init(pstm, filename, ispin)
    type(ps_tm_t),    intent(inout) :: pstm
    character(len=*), intent(in)    :: filename
    integer,          intent(in)    :: ispin

    character(len=256) :: filename2
    integer :: iunit, l, n
    logical :: found
    FLOAT :: x

    call push_sub('ps_tm.ps_tm_init')

    ! Sets the spin components
    pstm%ispin = ispin

    ! Find out where the hell the file is.
    filename2 = trim(filename) // '.vps'
    inquire(file=filename2, exist=found)
    message(1) = "Reading pseudopotential from file:"
    if(found) then
      write(message(2), '(6x,3a)') "'", trim(filename2), "'"
      call write_info(2)

      iunit = io_open(filename2, action='read', form='unformatted', status='old')
      call ps_tm_file_read(iunit, .false., pstm%tm_file)
      call io_close(iunit)
    else
      filename2 = trim(filename) // '.ascii'
      inquire(file=filename2, exist=found)
      if(.not.found) then
        filename2 = trim(conf%share) // "/PP/TM2/" // trim(filename) // ".ascii"
        inquire(file=filename2, exist=found)
        if(.not.found) then
          message(1) = "Pseudopotential file '" // trim(filename) // "{.vps|.ascii}' not found"
          call write_fatal(1)
        end if
      end if

      write(message(2), '(6x,3a)') "'", trim(filename2), "'"
      call write_info(2)

      iunit = io_open(filename2, action='read', form='formatted', status='old')
      call ps_tm_file_read(iunit, .true., pstm%tm_file)
      call io_close(iunit)
    end if

    ! Fills the valence configuration data.
    call build_valconf(pstm%tm_file, ispin, pstm%conf)

    ! Hack
    if(mod(pstm%tm_file%nr, 2) == 0) pstm%tm_file%nr = pstm%tm_file%nr - 1

    call file_to_grid(pstm%tm_file, pstm%ps_grid)

    call pop_sub()
  end subroutine ps_tm_init

  
  ! ---------------------------------------------------------
  subroutine ps_tm_end(ps_tm)
    type(ps_tm_t), intent(inout) :: ps_tm

    call ps_in_grid_end(ps_tm%ps_grid)
    call ps_tm_file_end(ps_tm%tm_file)
  end subroutine ps_tm_end


  ! ---------------------------------------------------------
  subroutine build_valconf(tm_file, ispin, conf)
    type(ps_tm_file_t), intent(in)  :: tm_file
    integer,            intent(in)  :: ispin
    type(valconf_t),    intent(out) :: conf

    character(len=1)  :: char1(6), char2
    integer :: i, l
    FLOAT   :: x

    call push_sub('tm.build_valconf')

    call valconf_null(conf)
    conf%symbol = tm_file%namatm
    conf%p = tm_file%npotd
    write(char2,'(i1)') conf%p

    select case(tm_file%irel)
    case('nrl')
      ! It seems that with some compiler version in the SP4 we can not have string
      ! concat inside the read statement, so we define read_fmt outside
      read(tm_file%title, '('//char2//'(i1,a1,f5.2,10x))')        &
        (conf%n(l), char1(l), conf%occ(l,1), l = 1, conf%p)
    case('isp')
      read(tm_file%title, '('//char2//'(i1,a1,f4.2,1x,f4.2,6x))') &
        (conf%n(l), char1(l), conf%occ(l,1), conf%occ(l,2), l = 1, conf%p)
    case('rel')
      read(tm_file%title, '('//char2//'(i1,a1,f5.2,10x))')        &
        (conf%n(l), char1(l), conf%occ(l,1), l = 1, conf%p)
    end select

    do l = 1, conf%p
      select case(char1(l))
      case('s'); conf%l(l) = 0
      case('p'); conf%l(l) = 1
      case('d'); conf%l(l) = 2
      case('f'); conf%l(l) = 3
      case default
        message(1) = 'Error reading pseudopotential file.'
        call write_fatal(1)
      end select
    end do

    do i = 1, conf%p
      l = conf%l(i)
      if(ispin==2 .and. tm_file%irel.ne.'isp') then
        x = conf%occ(i, 1)
        conf%occ(i, 1) = min(x, real(2*l+1, PRECISION))
        conf%occ(i, 2) = x - conf%occ(i,1)
      end if
    end do

    call pop_sub()
  end subroutine build_valconf 

  !----------------------------------------------------------------
  subroutine file_to_grid(tm_file, ps_grid)
    type(ps_tm_file_t), intent(in)  :: tm_file
    type(ps_in_grid_t), intent(out) :: ps_grid

    ! Initializes the pseudo in the logaritmic grid.
    call ps_in_grid_init(ps_grid,                      &
      LOGRID_TM, tm_file%a, tm_file%b, tm_file%nr,  &
      tm_file%npotd, tm_file%npotu)
    
    ps_grid%zval      = tm_file%zval
    ps_grid%vps(:,:)  = tm_file%vps(:,:)
    ps_grid%chcore(:) = tm_file%chcore(:)
    if(ps_grid%so_no_l_channels > 0) then
      ps_grid%so_vps(:,:) = tm_file%vso(:,:)
    end if

    ps_grid%core_corrections = .true.
    if(trim(tm_file%icore) == 'nc') ps_grid%core_corrections = .false.

  end subroutine file_to_grid


  ! ---------------------------------------------------------
  subroutine ps_tm_process(ps_tm, lmax, lloc)
    type(ps_tm_t), intent(inout) :: ps_tm
    integer,       intent(in)    :: lmax, lloc

    call push_sub('tm.tm_process')

    ! get the pseudoatomic eigenfunctions
    ALLOCATE(ps_tm%eigen(ps_tm%tm_file%npotd, 3), ps_tm%tm_file%npotd*3)
    call solve_schroedinger(ps_tm%tm_file, ps_tm%ps_grid%g, &
      ps_tm%conf, ps_tm%ispin, ps_tm%ps_grid%rphi, ps_tm%eigen)

    ! check norm of rphi
    call ps_in_grid_check_rphi(ps_tm%ps_grid)

    ! Fix the local potential. Final argument is the core radius
    call ps_in_grid_vlocal(ps_tm%ps_grid, lloc, M_THREE)

    ! Calculate kb cosines and norms
    call ps_in_grid_kb_cosines(ps_tm%ps_grid, lloc)

    ! Ghost analysis.
    call ghost_analysis(ps_tm%tm_file, ps_tm%ps_grid, ps_tm%ps_grid%g, ps_tm%eigen, lmax)

    ! Define the KB-projector cut-off radii
    call ps_in_grid_cutoff_radii(ps_tm%ps_grid, lloc)

    ! Calculate KB-projectors
    call ps_in_grid_kb_projectors(ps_tm%ps_grid)

    deallocate(ps_tm%eigen)
    call pop_sub()
  end subroutine ps_tm_process


  ! ---------------------------------------------------------
  subroutine solve_schroedinger(tm_file, g, conf, ispin, rphi, eigen)
    type(ps_tm_file_t), intent(inout) :: tm_file ! WARNING: should be intent(in)
    type(logrid_t),     intent(in)    :: g
    type(valconf_t),    intent(in)    :: conf
    integer,            intent(in)    :: ispin
    FLOAT,              intent(out)   :: rphi(:,:,:), eigen(:,:)
    

    character(len=3) :: functl
    integer :: iter, ir, is, l, nnode, nprin, ierr !, irel
    FLOAT :: vtot, diff, a2b4
    FLOAT, allocatable :: ve(:, :), rho(:, :), prev(:, :)
    FLOAT, parameter :: tol = CNST(1.0e-10)

    ! These variables are in double precision, no matter if single precision version of
    ! octopus is compiled, because they are passed to egofv.
    DOUBLE :: e, z, dr, rmax
    DOUBLE, allocatable :: s(:), hato(:), gg(:), y(:)

    call push_sub('tm.solve_schroedinger')

    ! Let us be a bit informative.
    message(1) = '      Calculating atomic pseudo-eigenfunctions for specie ' // tm_file%namatm // '....'
    call write_info(1)

    ! Allocation.
    ALLOCATE(s   (g%nrval),            g%nrval)
    ALLOCATE(hato(g%nrval),            g%nrval)
    ALLOCATE(gg  (g%nrval),            g%nrval)
    ALLOCATE(y   (g%nrval),            g%nrval)
    ALLOCATE(ve  (g%nrval, ispin), g%nrval*ispin)
    ALLOCATE(rho (g%nrval, ispin), g%nrval*ispin)
    ALLOCATE(prev(g%nrval, ispin), g%nrval*ispin)
    s = M_ZERO; hato = M_ZERO; gg = M_ZERO; y = M_ZERO;
    ve = M_ZERO; rho = M_ZERO; prev = M_ZERO

    ! These numerical parameters have to be fixed for egofv to work.
    s(:) = g%drdi(:)**2
    a2b4 = M_FOURTH*tm_file%a**2

    !  ionic pseudopotential if core-correction for hartree
    if((tm_file%icore == 'pche') .or. (tm_file%icore == 'fche')) then
      call vhrtre(tm_file%chcore, ve(:, 1), g%rofi, g%drdi, g%s, g%nrval, g%a)
      do l = 1, tm_file%npotd
        tm_file%vps(:, l) = tm_file%vps(:, l) + ve(:, 1)
      end do
    end if

    ! Calculation of the valence screening potential from the density:
    !       ve(1:nrval) is the hartree+xc potential created by the pseudo -
    !               valence charge distribution (everything in Rydberts, and bohrs)
    rho(:, 1) = tm_file%rho_val(:)
    select case(tm_file%icorr)
    case('pb')
      functl = 'GGA'
    case default
      functl = 'LDA'
    end select

    call atomhxc(functl, g, 1, rho(:, 1:1), ve(:, 1:1), tm_file%chcore)

    ! Calculation of the pseudo-wave functions.
    !       rphi(1:nrval, 1:conf%p, :) : radial pseudo-wave functions. They are normalized so that
    !           int dr {rphi^2} = 1. Thus its units are of bohr^(-1/2).
    !       eigen(1:conf%p, :)        : eigenvalues, in Rydbergs.
    do l = 1, tm_file%npotd
      do ir = 2, g%nrval
        vtot = tm_file%vps(ir, l) + ve(ir, 1) + dble((l-1)*l)/(g%rofi(ir)**2)
        hato(ir) = vtot*s(ir) + a2b4
      end do
      hato(1) = linear_extrapolate(g%rofi(1), g%rofi(2), g%rofi(3), &
        hato(2), hato(3))

      nnode = 1
      nprin = l
      e     = -((tm_file%zval/dble(nprin))**2)
      z     = tm_file%zval
      dr    = CNST(-1.0e5)
      rmax = g%rofi(g%nrval)

      call egofv(hato, s, g%nrval, e, gg, y, l, z, &
        real(g%a, 8), real(g%b, 8), rmax, nprin, nnode, dr, ierr)
      
      if(ierr.ne.0) then
        write(message(1),'(a)') 'The algorithm that calculates atomic wave functions could not'
        write(message(2),'(a)') 'do its job. The program will terminate, since the wavefunctions'
        write(message(3),'(a)') 'are needed. Change the pseudopotential or improve the code.'
        call write_fatal(3)
      end if
      eigen(l, 1) = e

      rphi(:, l, 1) = gg(:) * sqrt(g%drdi(:))

    end do

    ! Make a copy of the wavefunctions.
    rphi(:, :, 2) = rphi(:, :, 1)

    ! And now, for the spin polarized case...
    spin_polarized: if(ispin == 2) then

      rho = M_ZERO    ! Here information of previous calculation could be used, but
      prev = M_ZERO   ! to save code lines, let us start from scratch.
      diff = CNST(1.0e5)
      iter = 0
      self_consistent: do
        prev = rho
        iter = iter + 1

        spin: do is = 1, ispin
          ang: do l = 1, tm_file%npotd
            do ir = 2, g%nrval
              vtot = tm_file%vps(ir, l) + ve(ir, is) + dble((l-1)*l)/(g%rofi(ir)**2)
              hato(ir) = vtot*s(ir) + a2b4
            end do
            hato(1) = linear_extrapolate(g%rofi(1), g%rofi(2), g%rofi(3), &
              hato(2), hato(3))

            nnode = 1
            nprin = l
            e     = -((tm_file%zval/dble(nprin))**2)
            z     = tm_file%zval
            dr    = -CNST(1.0e5)
            rmax = g%rofi(g%nrval)

            call egofv(hato, s, g%nrval, e, gg, y, l, z, &
              real(g%a, 8), real(g%b, 8), rmax, nprin, nnode, dr, ierr)

            if(ierr.ne.0) then
              write(message(1),'(a)') 'The algorithm that calculates atomic wave functions could not'
              write(message(2),'(a)') 'do its job. The program will terminate, since the wavefunctions'
              write(message(3),'(a)') 'are needed. Change the pseudopotential or improve the code.'
              call write_fatal(3)
            end if
            eigen(l, 1 + is) = e

            rphi(:, l, 1 + is) = gg(:) * sqrt(g%drdi(:))
          end do ang
        end do spin

        rho = M_ZERO
        do is = 1, ispin
          do l = 1, tm_file%npotd
            rho(:, is) = rho(:, is) + conf%occ(l, is)*rphi(:, l, 1 + is)**2
          end do
        end do

        diff = M_ZERO
        do is = 1, ispin
          diff = diff + sqrt( sum(g%drdi(:) * (rho(:, is) - prev(:, is))**2) )
        end do

        if(diff < tol) exit self_consistent
        if(iter>1) rho = 0.5*rho + 0.5*prev

        write(message(1),'(a,i4,a,e10.2)') '      Iter =',iter,'; Diff =',diff
        call write_info(1)

        call atomhxc(functl, g, 2, rho, ve, tm_file%chcore)
      end do self_consistent

    end if spin_polarized

    ! Exit this...
    message(1) = '      Done.'; call write_info(1)
    deallocate(s, ve, hato, gg, y)

    call pop_sub()
  end subroutine solve_schroedinger


  ! ---------------------------------------------------------
  subroutine ghost_analysis(tm_file, ps_grid, g, eigen, lmax)
    type(ps_tm_file_t), intent(in) :: tm_file
    type(ps_in_grid_t), intent(in) :: ps_grid
    type(logrid_t),     intent(in) :: g
    FLOAT,              intent(in) :: eigen(:,:)
    integer,            intent(in) :: lmax

    character(len=3) :: functl
    integer :: ir, l, nnode, nprin, ighost, ierr
    FLOAT :: vtot, a2b4
    FLOAT, allocatable :: ve(:), elocal(:,:)

    DOUBLE :: z, e, dr, rmax
    DOUBLE, allocatable :: hato(:), s(:), gg(:), y(:)

    call push_sub('tm.ghost_analysis')

    ALLOCATE(ve    (g%nrval), g%nrval)
    ALLOCATE(s     (g%nrval), g%nrval)
    ALLOCATE(hato  (g%nrval), g%nrval)
    ALLOCATE(gg    (g%nrval), g%nrval)
    ALLOCATE(y     (g%nrval), g%nrval)
    ALLOCATE(elocal(2, lmax+1),  2*(lmax+1))

    select case(tm_file%icorr)
    case('pb')
      functl = 'GGA'
    case default
      functl = 'LDA'
    end select

    call atomhxc(functl, g, 1, tm_file%rho_val(:), ve(:), ps_grid%chcore(:))

    ! calculate eigenvalues of the local potential for ghost analysis
    s(:) = g%drdi(:)**2
    a2b4 = M_FOURTH*g%a**2

    do l = 1, lmax+1
      do ir = 2, g%nrval
        vtot = ps_grid%vlocal(ir) + ve(ir) + dble((l-1)*l)/(g%rofi(ir)**2)
        hato(ir) = vtot*s(ir) + a2b4
      end do
      hato(1) = linear_extrapolate(g%rofi(1), g%rofi(2), g%rofi(3), &
        hato(2), hato(3))      

      do nnode = 1, 2
        nprin = l
        e     = -(ps_grid%zval/dble(nprin))**2
        z     = ps_grid%zval
        dr    = CNST(-1.0e5)
        rmax  = g%rofi(g%nrval)

        call egofv(hato, s, g%nrval, e, gg, y, l, z, &
          real(g%a, 8), real(g%b, 8), rmax, nprin, nnode, dr, ierr)

        elocal(nnode, l) = e
      end do
    end do

    ! Ghost analysis
    do l = 1, lmax+1
      ighost = -1

      if(ps_grid%dkbcos(l) > M_ZERO) then
        if(eigen(l, 1) > elocal(2, l)) then
          ighost = 1
        end if
      else if(ps_grid%dkbcos(l) < M_ZERO) then
        if(eigen(l, 1) > elocal(1, l)) then
          ighost = 1
        end if
      end if

      if(ighost >= 0) then
        write(message(1), '(a,i2)') "Ghost state found for l = ", l
        call write_warning(1)
      end if
    end do

    deallocate(hato, gg, y, elocal, s, ve)

    call pop_sub()
  end subroutine ghost_analysis

end module ps_tm_m
