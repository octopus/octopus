!! Copyright (C) 2002-2011 M. Marques, A. Castro, A. Rubio, G. Bertsch, M. Oliveira
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

module ps_upf_m
  use atomic_m
  use global_m
  use io_m
  use messages_m
  use profiling_m
  use ps_in_grid_m

  implicit none

  private
  public ::     &
    ps_upf_t,       &
    ps_upf_init,    &
    ps_upf_end

  type ps_upf_t
    type(valconf_t)    :: conf

    integer :: kb_nc
    integer :: l_local
    FLOAT :: local_radius
    FLOAT, pointer :: kb_radius(:)

    ! The contents of the file:
    !Header
    integer      :: version  ! UPF file version number
    character(3) :: symbol   ! Element label
    character(2) :: type     ! Pseudo type (NC or US). In octopus only NC is implemented.
    FLOAT        :: z_val    ! z valence
    integer      :: l_max    ! maximum angular momentum component
    integer      :: n_proj   ! number of projectors
    integer      :: n_wfs    ! number of wavefunctions
    integer, pointer :: n(:)
    integer, pointer :: l(:)
    FLOAT, pointer :: occ(:)

    !Radial mesh
    integer        :: np
    FLOAT, pointer :: r(:)
    FLOAT, pointer :: drdi(:)

    !nlcc
    logical        :: nlcc
    FLOAT, pointer :: core_density(:)

    !KB projectors
    FLOAT,   pointer :: v_local(:)
    integer, pointer :: proj_l(:)
    integer, pointer :: proj_np(:)
    FLOAT,   pointer :: proj(:,:)
    FLOAT,   pointer :: e(:)

    !Wavefunctions
    FLOAT,  pointer :: wfs(:,:)

    !Valence charge
    FLOAT,  pointer :: rho(:)
    
    !Extra information
    FLOAT, pointer :: proj_j(:)

  end type ps_upf_t

contains

  ! ---------------------------------------------------------
  subroutine ps_upf_init(ps_upf, filename)
    type(ps_upf_t),   intent(inout) :: ps_upf
    character(len=*), intent(in)    :: filename

    character(len=256) :: filename2
    integer :: iunit, l
    logical :: found
    logical, allocatable :: found_l(:)    

    PUSH_SUB(ps_upf_init)

    ! Find out the file and read it.
    filename2 = trim(filename) // '.UPF'
    inquire(file=filename2, exist=found)
    message(1) = "Reading pseudopotential from file:"

    if(.not. found) then
      filename2 = trim(conf%share) // "/PP/UPF/" // trim(filename) // ".UPF"
      inquire(file=filename2, exist=found)

      if(.not.found) then
        message(1) = "Pseudopotential file '" // trim(filename) // ".UPF' not found"
        call messages_fatal(1)
      end if
    end if

    write(message(2), '(6x,3a)') "'", trim(filename2), "'"
    call messages_info(2)
    iunit = io_open(filename2, action='read', form='formatted', status='old', is_tmp=.true.)
    call ps_upf_file_read(iunit, ps_upf)
    call io_close(iunit)

    !Build valence configuration
    call valconf_null(ps_upf%conf)
    ps_upf%conf%symbol = ps_upf%symbol
    ps_upf%conf%type = 0
    ps_upf%conf%p = ps_upf%n_wfs
    ps_upf%conf%n(1:ps_upf%n_wfs) = ps_upf%n(1:ps_upf%n_wfs)
    ps_upf%conf%l(1:ps_upf%n_wfs) = ps_upf%l(1:ps_upf%n_wfs)
    ps_upf%conf%occ(1:ps_upf%n_wfs,1) = ps_upf%occ(1:ps_upf%n_wfs)

    ps_upf%l_max = maxval(ps_upf%l)

    !Check if the local component is one of the angular momentum channels
    SAFE_ALLOCATE(found_l(0:ps_upf%l_max))
    found_l = .true.
    do l = 0, ps_upf%l_max
      if (any(ps_upf%proj_l == l)) then
        found_l(l) = .false.
      end if
    end do
    if (count(found_l) /= 1) then
      ps_upf%l_local = -1
    else
      do l = 0, ps_upf%l_max
        if (found_l(l)) then
          ps_upf%l_local = l
          exit
        end if
      end do
    end if

    write(message(1), '(a,i2)') '      l max = ', ps_upf%l_max
    write(message(2), '(a,i2)') '      l loc = ', ps_upf%l_local
    call messages_info(2)

    ! Define the KB-projector cut-off radii
    call ps_upf_cutoff_radii(ps_upf)

    ! check norm of rphi
    call ps_upf_check_rphi(ps_upf)

    POP_SUB(ps_upf_init)
  end subroutine ps_upf_init

  ! ---------------------------------------------------------
  subroutine ps_upf_end(ps_upf)
    type(ps_upf_t), intent(inout) :: ps_upf

    PUSH_SUB(ps_upf_end)

    SAFE_DEALLOCATE_P(ps_upf%kb_radius)

    SAFE_DEALLOCATE_P(ps_upf%n)
    SAFE_DEALLOCATE_P(ps_upf%l)
    SAFE_DEALLOCATE_P(ps_upf%occ)
    SAFE_DEALLOCATE_P(ps_upf%core_density)
    SAFE_DEALLOCATE_P(ps_upf%r)
    SAFE_DEALLOCATE_P(ps_upf%drdi)
    SAFE_DEALLOCATE_P(ps_upf%v_local)
    SAFE_DEALLOCATE_P(ps_upf%proj_l)
    SAFE_DEALLOCATE_P(ps_upf%proj_np)
    SAFE_DEALLOCATE_P(ps_upf%proj)
    SAFE_DEALLOCATE_P(ps_upf%e)
    SAFE_DEALLOCATE_P(ps_upf%wfs)
    SAFE_DEALLOCATE_P(ps_upf%rho)
    SAFE_DEALLOCATE_P(ps_upf%proj_j)

    POP_SUB(ps_upf_end)
  end subroutine ps_upf_end

  ! ---------------------------------------------------------
  subroutine ps_upf_file_read(unit, ps_upf)
    integer,        intent(in)    :: unit
    type(ps_upf_t), intent(inout) :: ps_upf

    integer :: ip, np, i, ir, idummy, ii, jj, n_dij
    character (len=2) :: nl
    character (len=80) :: dummy
    logical :: ok
    FLOAT :: fdummy

    PUSH_SUB(ps_upf_file_read)

    ps_upf%kb_nc = 1

    !Header info
    call init_tag(unit, "PP_HEADER", .true.)

    read(unit,*) ps_upf%version, dummy  ! n        "Version Number"
    read(unit,*) ps_upf%symbol, dummy   ! psd      "Element"
    read(unit,*) ps_upf%type, dummy     ! US|NC|PAW    "Ultrasoft|Norm conserving|Projector-augmented"
    if (ps_upf%type /= "NC") then
      message(1) = "Octopus can only read norm-conserving pseudo-potentials from UPF format."
      call messages_fatal(1)
    end if
    read(unit,*) ps_upf%nlcc, dummy   ! nlcc     "Nonlinear Core Correction"
    read(unit,*) dummy                ! dft      "Exch-Corr"
    read(unit,*) ps_upf%z_val, dummy  ! zp       "Z valence"
    read(unit,*) dummy                ! etotps   "Total Energy"
    read(unit,*) dummy                ! ecutwfc, ecutrho     "Suggested Cutoff for wfc and rho"
    read(unit,*) dummy                ! lmax     "Max angular momentum component", THIS IS NOT THE LMAX WE NEED
    read(unit,*) np, dummy            ! mesh     "Number of points in mesh"
    read(unit,*) ps_upf%n_wfs, ps_upf%n_proj, dummy !  natwfc, nbeta   "Number of wavefunctions, projectors"
    read(unit,*) dummy                ! "Wavefunctions   nl   l   occ"
    SAFE_ALLOCATE(ps_upf%n(1:ps_upf%n_wfs))
    SAFE_ALLOCATE(ps_upf%l(1:ps_upf%n_wfs))
    SAFE_ALLOCATE(ps_upf%occ(1:ps_upf%n_wfs))

    ! els(1)      lchi(1)      oc(1)
    !  ...
    ! els(natwfc) lchi(natwfc) oc(natwfc)
    do i = 1, ps_upf%n_wfs
      read(unit,*) nl, ps_upf%l(i), ps_upf%occ(i)
      ps_upf%n(i) = iachar(nl(1:1)) - iachar('0')
      ! some pseudopotentials do not have a number, but just a letter,
      ! so we assume the level is 1a
      if(ps_upf%n(i) < 1 .or. ps_upf%n(i) > 9)  ps_upf%n(i) = 1
    end do

    call check_end_tag(unit, "PP_HEADER")

    !Mesh info
    call init_tag(unit, "PP_MESH", .true.)

    call init_tag(unit, "PP_R", .false.)
    read(unit,*) fdummy
    if (fdummy /= M_ZERO) then
      ps_upf%np = np + 1
      ip = 2
    else
      ps_upf%np = np
      ip = 1
    end if
    SAFE_ALLOCATE(ps_upf%r(1:ps_upf%np))
    ps_upf%r(1) = M_ZERO
    ps_upf%r(ip) = fdummy
    call init_tag(unit, "PP_R", .true.)
    read(unit,*) (ps_upf%r(ir), ir = ip, ps_upf%np)
    call check_end_tag(unit, "PP_R")

    call init_tag(unit, "PP_RAB", .false.)
    SAFE_ALLOCATE(ps_upf%drdi(1:ps_upf%np))
    read(unit,*) (ps_upf%drdi(ir), ir = ip, ps_upf%np)
    call check_end_tag(unit, "PP_RAB")
    if (ip == 2) ps_upf%drdi(1) = M_ZERO
    call check_end_tag(unit, "PP_MESH")

    !Non-linear core-corrections
    if (ps_upf%nlcc) then
      call init_tag(unit, "PP_NLCC", .true.)
      SAFE_ALLOCATE(ps_upf%core_density(1:ps_upf%np))
      read(unit,*) (ps_upf%core_density(ir), ir = ip, ps_upf%np)
      call check_end_tag(unit, "PP_NLCC")
      if (ip == 2) ps_upf%core_density(1) = linear_extrapolate(ps_upf%r(1), ps_upf%r(2), &
           ps_upf%r(3), ps_upf%core_density(2), ps_upf%core_density(3))
    else
      nullify(ps_upf%core_density)
    end if

    !Local component
    call init_tag(unit, "PP_LOCAL", .true.)
    SAFE_ALLOCATE(ps_upf%v_local(1:ps_upf%np))
    read(unit,*) (ps_upf%v_local(ir), ir = ip, ps_upf%np)
    if (ip == 2) ps_upf%v_local(1) = linear_extrapolate(ps_upf%r(1), ps_upf%r(2), &
           ps_upf%r(3), ps_upf%v_local(2), ps_upf%v_local(3))
    call check_end_tag(unit, "PP_LOCAL")

    !Non-local components
    call init_tag(unit, "PP_NONLOCAL", .true.)

    SAFE_ALLOCATE(ps_upf%proj(1:ps_upf%np, 1:ps_upf%n_proj))
    SAFE_ALLOCATE(ps_upf%proj_l(1:ps_upf%n_proj))
    SAFE_ALLOCATE(ps_upf%proj_np(1:ps_upf%n_proj))
    ps_upf%proj = M_ZERO
    do i = 1, ps_upf%n_proj
      call init_tag(unit, "PP_BETA", .false.)
      read(unit,*) idummy, ps_upf%proj_l(i), dummy
      read(unit,*) ps_upf%proj_np(i)
      read(unit,*) (ps_upf%proj(ir, i), ir = ip, ps_upf%proj_np(i)+ip-1)
      if (ip == 2) ps_upf%proj(1, i) = M_ZERO!linear_extrapolate(ps_upf%r(1), ps_upf%r(2), &
      !ps_upf%r(3), ps_upf%proj(2, i), ps_upf%proj(3, i))
      call check_end_tag(unit, "PP_BETA", ok = ok)
      if(.not. ok) then
        ! This is a 'new' UPF file generated by Quantum Espresso that
        ! is not compatible with the Quantum Espresso
        ! specifications. We have to skip one line (one line was read
        ! by check_end_tag already).
        read(unit,*) dummy
        call check_end_tag(unit, "PP_BETA")
      end if
    end do
    
    call init_tag(unit, "PP_DIJ", .false.)
    SAFE_ALLOCATE(ps_upf%e(1:ps_upf%n_proj))
    ps_upf%e = M_ZERO
    read(unit,*) n_dij, dummy
    do i = 1, n_dij
      read(unit,*) ii, jj, ps_upf%e(ii)
      if (ii /= jj) then
        message(1) = "Error while reading pseudo-potential data."
        call messages_fatal(1)
      end if
    end do
    call check_end_tag(unit, "PP_DIJ")

    call check_end_tag(unit, "PP_NONLOCAL")

    !Pseudo wavefunctions
    call init_tag(unit, "PP_PSWFC", .true.)
    SAFE_ALLOCATE(ps_upf%wfs(1:ps_upf%np, 1:ps_upf%n_wfs))
    do i = 1, ps_upf%n_wfs
      read(unit,*) dummy
      read(unit,*) (ps_upf%wfs(ir, i), ir = ip, ps_upf%np)
      if (ip == 2) ps_upf%wfs(1, i) = M_ZERO!linear_extrapolate(ps_upf%r(1), ps_upf%r(2), &
      !ps_upf%r(3), ps_upf%wfs(2, i), ps_upf%wfs(3, i))
    end do
    call check_end_tag(unit, "PP_PSWFC")

    !Valence charge
    call init_tag(unit, "PP_RHOATOM", .true.)
    SAFE_ALLOCATE(ps_upf%rho(1:ps_upf%np))
    read(unit,*) (ps_upf%rho(ir), ir = ip, ps_upf%np)
      if (ip == 2) ps_upf%rho(1) = M_ZERO!linear_extrapolate(ps_upf%r(1), ps_upf%r(2), &
      !ps_upf%r(3), ps_upf%rho(2), ps_upf%rho(3))
    call check_end_tag(unit, "PP_RHOATOM")

    !Extra information
    if (tag_isdef(unit, "PP_ADDINFO")) then
      ps_upf%kb_nc = 2
      call init_tag(unit, "PP_ADDINFO", .true.)
      do i = 1, ps_upf%n_wfs
        read(unit,*) dummy
      end do
      SAFE_ALLOCATE(ps_upf%proj_j(1:ps_upf%n_proj))
      do i = 1, ps_upf%n_proj
        read(unit,*) dummy, ps_upf%proj_j(i)
      end do
      read(unit,*) dummy
      call check_end_tag(unit, "PP_ADDINFO")
    else
      nullify(ps_upf%proj_j)
    end if

    POP_SUB(ps_upf_file_read)
  end subroutine ps_upf_file_read

  subroutine init_tag(unit, string, go_back)
    integer,          intent(in) :: unit
    character(len=*), intent(in) :: string
    logical,          intent(in) :: go_back

    integer :: iostat
    character(len=80) :: string2

    PUSH_SUB(init_tag)

    if (go_back) rewind(unit)
    do
      read(unit, *, iostat = iostat) string2
      if (iostat > 0) then
        message(1) = "Error in subroutine init_tag"
        call messages_fatal(1)
      elseif (iostat == -1) then
        write(message(1),*) "No ", trim(string), " tag found."
        call messages_fatal(1)
      end if
      if (string_matches("<"//string//">", string2) ) exit
    end do
write(*,*) string    
    POP_SUB(init_tag)
  end subroutine init_tag

  logical function tag_isdef(unit, string)
    integer,          intent(in) :: unit
    character(len=*), intent(in) :: string
      
    integer :: iostat
    character(len=80) :: string2
    
    PUSH_SUB(tag_isdef)

    rewind(unit)
    do
      read(unit, *, iostat = iostat) string2
      if (iostat > 0) then
        message(1) = "Error in subroutine tag_isdef"
        call messages_fatal(1)
      elseif (iostat < 0) then
        tag_isdef = .false.
        exit
      end if
      if (string_matches("<"//string//">", string2) ) then
        tag_isdef = .true.
        exit
      end if
    end do

    POP_SUB(tag_isdef)
  end function tag_isdef
  
  ! ------------------------------------------------------------------------
  
  subroutine check_end_tag(unit, string, ok)
    integer,           intent(in)  :: unit
    character(len=*),  intent(in)  :: string
    logical, optional, intent(out) :: ok

    integer :: iostat
    character(len=80) :: string2
      
    PUSH_SUB(check_end_tag)

    read(unit, '(a)', iostat = iostat) string2
    if ((.not. string_matches("</"//string//">", string2)) .or. iostat /= 0) then
      if(present(ok)) then
        ok = .false.
      else
        write(message(1),*) "Could not find closing tag </", string, ">."
        call messages_fatal(1)
      end if
    else
      if(present(ok)) ok = .true.
    end if

    POP_SUB(check_end_tag)
  end subroutine check_end_tag

  ! ------------------------------------------------------------------------

  logical function string_matches(string1, string2)
    character(len=*), intent(in) :: string1, string2
      
    integer :: l1, l2, l
    
    string_matches = .false.

    l1 = len_trim(string1)
    l2 = len_trim(string2)
    do l = 1, (l2 - l1 + 1)
      if (string1(1:l1) == string2(l:(l+l1-1))) then
        string_matches = .true.
        exit
      end if
    end do
      
  end function string_matches

  ! ---------------------------------------------------------
  subroutine ps_upf_cutoff_radii(ps_upf)
    type(ps_upf_t), intent(inout) :: ps_upf

    integer          :: i, ir
    FLOAT            :: dincv
    FLOAT, parameter :: threshold = CNST(1.0e-6)

    PUSH_SUB(ps_upf_cutoff_radii)

    ! local part ....
    do ir = ps_upf%np-1, 2, -1
      dincv = abs(ps_upf%v_local(ir)*ps_upf%r(ir) + M_TWO*ps_upf%z_val)
      if(dincv > threshold) exit
    end do
    ps_upf%local_radius = ps_upf%r(ir + 1)

    ! non-local part....
    SAFE_ALLOCATE(ps_upf%kb_radius(1:ps_upf%n_proj))

    do i = 1, ps_upf%n_proj

      do ir = ps_upf%np-1, 2, -1
        dincv = abs(ps_upf%proj(ir, i))

        if(dincv > threshold) exit
      end do

      ps_upf%kb_radius(i) = ps_upf%r(ir + 1)
    end do
    
    POP_SUB(ps_upf_cutoff_radii)
  end subroutine ps_upf_cutoff_radii


  ! ---------------------------------------------------------
  ! checks normalization of the pseudo wavefunctions
  subroutine ps_upf_check_rphi(ps_upf)
    type(ps_upf_t), intent(in) :: ps_upf
    
    integer :: i
    FLOAT   :: nrm

    PUSH_SUB(ps_upf_check_rphi)

    !  checking normalization of the wavefunctions
    do i = 1, ps_upf%n_wfs
      nrm = sqrt(sum(ps_upf%drdi(:)*ps_upf%wfs(:, i)**2))
      nrm = abs(nrm - M_ONE)
      if (nrm > CNST(1.0e-5)) then
        write(message(1), '(a,i2,a)') "Eigenstate for l = ", ps_upf%l(i), ' is not normalized'
        write(message(2), '(a, f12.6,a)') '(abs(1 - norm) = ', nrm, ')'
        call messages_warning(2)
      end if
    end do

    POP_SUB(ps_upf_check_rphi)
  end subroutine ps_upf_check_rphi

end module ps_upf_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
