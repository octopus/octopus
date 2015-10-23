!! Copyright (C) 2014 M. Marques, A. Castro, A. Rubio, G. Bertsch, J.Jornet-Somoza
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
!! $Id: local_write.F90 11872 2014-03-12 19:43:35Z dstrubbe $

#include "global.h"

module local_write_m
  use box_union_m
  use iso_c_binding
  use geometry_m
  use global_m
  use grid_m
  use hamiltonian_m
  use io_m
  use io_function_m
  use kick_m
  use mesh_function_m
  use mesh_m
  use messages_m
  use mpi_m
  use mpi_debug_m
  use mpi_lib_m
  use multicomm_m
  use parser_m
  use poisson_m
  use profiling_m
  use species_m
  use states_m
  use unit_m
  use unit_system_m
  use varinfo_m
  use v_ks_m
  use write_iter_m

  implicit none

  private
  public ::         &
    local_write_t,     &
    local_write_check_hm, &
    local_write_init,  &
    local_write_iter,  &
    local_write_end   

  integer, parameter ::   &
    LOCAL_OUT_MULTIPOLES  =  1, &
    LOCAL_OUT_DENSITY     =  2, &
    LOCAL_OUT_POTENTIAL   =  3, &
    LOCAL_OUT_ENERGY      =  8, &
    LOCAL_OUT_MAX         =  8
  
  type local_write_prop_t
    private
    type(c_ptr) :: handle
    logical :: write = .false.
  end type local_write_prop_t

  type local_write_t
    private
    type(local_write_prop_t),allocatable :: out(:,:)
    integer                  :: how              !< how to output
    integer                  :: lmax             !< maximum multipole moment to output
  end type local_write_t

contains

  ! ---------------------------------------------------------
  logical function local_write_check_hm(writ) result(check)
    type(local_write_t), intent(out)   :: writ

    PUSH_SUB(local_write_check_hm)
    check = .false.
    if (any(writ%out(LOCAL_OUT_POTENTIAL, :)%write) .or. &
        any(writ%out(LOCAL_OUT_ENERGY, :)%write) ) check = .true.

    PUSH_SUB(local_write_check_hm)
  end function local_write_check_hm

  ! ---------------------------------------------------------
  subroutine local_write_init(writ, nd, lab, iter, dt)
    type(local_write_t), intent(out)   :: writ
    integer,             intent(in)    :: nd 
    character(len=15),   intent(in)    :: lab(:)
    integer,             intent(in)    :: iter
    FLOAT,               intent(in)    :: dt

    integer :: first, id, flags, iout, default

    PUSH_SUB(local_write_init)

    ! FIXME: if and when these routines are called from a normal run, the Section can be Output.
    ! but then it will need to be in a different folder, since src/util is not linked by the other folders.

    !%Variable LDOutput
    !%Type flag
    !%Default multipoles 
    !%Section Utilities::oct-local_multipoles
    !%Description
    !% Defines what should be output during the local domains 
    !% simulation. Many of the options can increase the computational
    !% cost of the simulation, so only use the ones that you need. In
    !% most cases the default value is enough, as it is adapted to the
    !% details. 
    !%Option multipoles 1
    !% Outputs the (electric) multipole moments of the density to the file <tt>ld.general/multipoles</tt>.
    !% This is required to, <i>e.g.</i>, calculate optical absorption spectra of finite systems. The
    !% maximum value of <math>l</math> can be set with the variable <tt>LDMultipoleLmax</tt>.
    !%Option density 2
    !% If set, <tt>octopus</tt> outputs the densities corresponding to the local domains to 
    !% the folder <tt>ld.general/densities</tt>. 
    !% The output format is set by the <tt>LDOutputFormat</tt> input variable.
    !%Option local_v 4
    !% If set, <tt>octopus</tt> outputs the different components of the potential
    !% to the folder <tt>ld.general/potential</tt>. 
    !% The output format is set by the <tt>LDOutputFormat</tt> input variable.
    !%Option energy 128
    !% If set, <tt>octopus</tt> outputs the different components of the energy of the local domains
    !% to the folder <tt>ld.general/energy</tt>.
    !%End

    default = 2**(LOCAL_OUT_MULTIPOLES - 1) 

    call parse_variable('LDOutput', default, flags)

    if(.not.varinfo_valid_option('LDOutput', flags, is_flag = .true.)) call messages_input_error('LDOutput')

    SAFE_ALLOCATE(writ%out(1:LOCAL_OUT_MAX, 1:nd))
    do iout = 1, LOCAL_OUT_MAX
      writ%out(iout,:)%write = (iand(flags, 2**(iout - 1)) /= 0)
    end do

    call messages_obsolete_variable('LDOutputHow', 'LDOutputFormat')
    
    !%Variable LDOutputFormat
    !%Type flag
    !%Default none
    !%Section Utilities::oct-local_multipoles
    !%Description
    !% Describes the format of the output files (see <tt>LDOutput</tt>).
    !% It can take the same values as <tt>OutputFormat</tt> flag.
    !%End
    call parse_variable('LDOutputFormat', 0, writ%how)
    if(.not.varinfo_valid_option('OutputFormat', writ%how, is_flag=.true.)) then
      call messages_input_error('LDOutputFormat')
    end if

    !%Variable LDMultipoleLmax
    !%Type integer
    !%Default 1
    !%Section Utilities::oct-local_multipoles
    !%Description
    !% Maximum electric multipole of the density output to the file <tt>local.multipoles/<>domain%<>.multipoles</tt>
    !% during a time-dependent simulation. Must be non-negative.
    !%End

    call parse_variable('LDMultipoleLmax', 1, writ%lmax)
    if (writ%lmax < 0) then
      write(message(1), '(a,i6,a)') "Input: '", writ%lmax, "' is not a valid LDMultipoleLmax."
      message(2) = '(Must be LDMultipoleLmax >= 0 )'
      call messages_fatal(2)
    end if

    call io_mkdir('local.general')

    if(mpi_grp_is_root(mpi_world)) then
      do id = 1, nd
        if(writ%out(LOCAL_OUT_MULTIPOLES, id)%write) then 
          call io_mkdir('local.general/multipoles')
          call write_iter_init(writ%out(LOCAL_OUT_MULTIPOLES,id)%handle, &
            iter, units_from_atomic(units_out%time, dt), &
          trim(io_workpath("local.general/multipoles/"//trim(lab(id))//".multipoles")))
        end if

        if(writ%out(LOCAL_OUT_POTENTIAL, id)%write) then
          call io_mkdir('local.general/potential')
          call write_iter_init(writ%out(LOCAL_OUT_POTENTIAL,id)%handle, first, &
          units_from_atomic(units_out%time, dt), trim(io_workpath("local.general/potential/"//trim(lab(id))//".potential")))
        end if

        if(writ%out(LOCAL_OUT_DENSITY, id)%write) then
          call io_mkdir('local.general/densities')
          call write_iter_init(writ%out(LOCAL_OUT_DENSITY,id)%handle, first, &
            units_from_atomic(units_out%time, dt), trim(io_workpath("local.general/densities/"//trim(lab(id))//".densities")))
        end if

        if(writ%out(LOCAL_OUT_ENERGY, id)%write) then 
          call io_mkdir('local.general/energy')
          call write_iter_init(writ%out(LOCAL_OUT_ENERGY,id)%handle, &
            iter, units_from_atomic(units_out%time, dt), &
          trim(io_workpath("local.general/energy/"//trim(lab(id))//".energy")))
        end if
      end do
    end if

    POP_SUB(local_write_init)
  end subroutine local_write_init

  ! ---------------------------------------------------------
  subroutine local_write_end(writ, nd)
    type(local_write_t), intent(inout) :: writ
    integer                            :: nd

    integer :: id, iout
    
    PUSH_SUB(local_write_end)

    do id = 1, nd
      do iout = 1, LOCAL_OUT_MAX
        if(writ%out(iout, id)%write)  call write_iter_end(writ%out(iout, id)%handle)
      end do
    end do
    SAFE_DEALLOCATE_A(writ%out)

    POP_SUB(local_write_end)
  end subroutine local_write_end

  ! ---------------------------------------------------------
  subroutine local_write_iter(writ, nd, domain, lab, inside, center, gr, st, & 
                              hm, ks, mc, geo, kick, iter, dt, l_start, l_end, ldoverwrite)
    type(local_write_t),    intent(inout) :: writ
    integer,                intent(in)    :: nd 
    type(box_union_t),      intent(in)    :: domain(:)
    character(len=15),      intent(in)    :: lab(:)
    logical,                intent(in)    :: inside(:,:)
    FLOAT  ,                intent(in)    :: center(:,:)
    type(grid_t),           intent(inout) :: gr
    type(states_t),         intent(inout) :: st
    type(hamiltonian_t),    intent(inout) :: hm
    type(v_ks_t),           intent(inout) :: ks
    type(multicomm_t),      intent(in)    :: mc
    type(geometry_t),       intent(inout) :: geo
    type(kick_t),           intent(inout) :: kick
    integer,                intent(in)    :: iter
    FLOAT,                  intent(in)    :: dt
    integer,                intent(in)    :: l_start
    integer,                intent(in)    :: l_end
    logical,                intent(in)    :: ldoverwrite

    type(profile_t), save :: prof
    integer :: id

    PUSH_SUB(local_write_iter)
    call profiling_in(prof, "LOCAL_WRITE_ITER")

    if(any(writ%out(LOCAL_OUT_MULTIPOLES,:)%write)) then
      if(.not.ldoverwrite)then
      end if
      call local_write_multipole(writ%out(LOCAL_OUT_MULTIPOLES, :), nd, domain, lab, inside, center, & 
                        gr, geo, st, writ%lmax, kick, iter, l_start, ldoverwrite, writ%how)
      do id = 1, nd
        call write_iter_flush(writ%out(LOCAL_OUT_MULTIPOLES, id)%handle)
      end do
    end if

    if(any(writ%out(LOCAL_OUT_DENSITY,:)%write).or.any(writ%out(LOCAL_OUT_POTENTIAL,:)%write)) &
      call local_write_density(writ%out(LOCAL_OUT_DENSITY, :), writ%out(LOCAL_OUT_POTENTIAL,:), & 
                               nd, domain, lab, inside, center, &
                               gr, geo, st, hm, ks, mc, iter, l_start, ldoverwrite, writ%how)
    
    if(any(writ%out(LOCAL_OUT_ENERGY, :)%write)) then
      call local_write_energy(writ%out(LOCAL_OUT_ENERGY, :), nd, domain, lab, inside, center, &
                               gr, geo, st, hm, ks, mc, iter, l_start, ldoverwrite, writ%how)
      do id = 1, nd
        call write_iter_flush(writ%out(LOCAL_OUT_ENERGY, id)%handle)
      end do
    end if

    call profiling_out(prof)
    POP_SUB(local_write_iter)
  end subroutine local_write_iter

  ! ---------------------------------------------------------
  subroutine local_write_density(out_dens, out_pot, nd, domain, lab, inside, center, & 
                                gr, geo, st, hm, ks, mc, iter, l_start, start, how)
    type(local_write_prop_t),      intent(inout) :: out_dens(:)
    type(local_write_prop_t),      intent(inout) :: out_pot(:)
    integer,                  intent(in)    :: nd 
    type(box_union_t),        intent(in)    :: domain(:)
    character(len=15),        intent(in)    :: lab(:)
    logical,                  intent(in)    :: inside(:,:)
    FLOAT,                    intent(in)    :: center(:,:)
    type(grid_t),         intent(inout) :: gr
    type(geometry_t),     intent(inout) :: geo
    type(states_t),       intent(inout) :: st
    type(hamiltonian_t),  intent(inout) :: hm
    type(v_ks_t),         intent(inout) :: ks
    type(multicomm_t),    intent(in)    :: mc
    integer,              intent(in) :: iter
    integer,              intent(in) :: l_start
    logical,              intent(in) :: start
    integer,              intent(in) :: how

    integer            :: id, is, ix, ierr
    character(len=120) :: folder, out_name
    FLOAT, allocatable :: tmp_rho(:),st_rho(:)
    FLOAT, allocatable :: tmp_vh(:)

    PUSH_SUB(local_write_density)

    SAFE_ALLOCATE(tmp_rho(1:gr%mesh%np))
    SAFE_ALLOCATE(st_rho(1:gr%mesh%np_part))
    SAFE_ALLOCATE(tmp_vh(1:gr%mesh%np))

    ! FIXME: use just v_ks_calc subroutine for either compute vhartree and vxc
    do id = 1, nd
      do is = 1, st%d%nspin
        tmp_rho = M_ZERO
        st_rho(:) = st%rho(:, is)
        do ix = 1, gr%mesh%np 
          if (inside(ix, id)) tmp_rho(ix) = st_rho(ix)
        end do
        if (out_dens(id)%write) then
          folder = 'local.general/densities/'//trim(lab(id))//'.densities/'
          write(out_name, '(a,a1,i0,a1,i7.7)')trim(lab(id)),'.',is,'.',iter
          call dio_function_output(how, &
            trim(folder), trim(out_name), gr%mesh, tmp_rho(1:gr%mesh%np), units_out%length, ierr, geo = geo)
        end if
        if (out_pot(id)%write) then
        !Computes Hartree potential just for n[r], r belongs to id domain.
          tmp_vh = M_ZERO
          call dpoisson_solve(psolver, tmp_vh, tmp_rho)
          folder = 'local.general/potential/'//trim(lab(id))//'.potential/'
          write(out_name, '(a,i0,a1,i7.7)')'vh.',is,'.',iter
          call dio_function_output(how, &
            trim(folder), trim(out_name), gr%mesh, tmp_vh, units_out%length, ierr, geo = geo)
        !Computes XC potential
          st%rho(:,is) = M_ZERO
          do ix = 1, gr%mesh%np
            if (inside(ix, id)) st%rho(ix, is) = st_rho(ix)
          end do
          call v_ks_calc(ks, hm, st, geo, calc_eigenval = .false. , calc_berry = .false. , calc_energy = .false.)
          folder = 'local.general/potential/'//trim(lab(id))//'.potential/'
          write(out_name, '(a,i0,a1,i7.7)')'vxc.',is,'.',iter
          call dio_function_output(how, &
            trim(folder), trim(out_name), gr%mesh, hm%vxc(:,is), units_out%length, ierr, geo = geo)
          st%rho(:,is) = st_rho(:)
        end if
      end do
    end do

    if (any(out_pot(:)%write)) then
      call v_ks_init(ks, gr, st, geo, mc)
      do is = 1, st%d%nspin
      !Computes Hartree potential
        call dpoisson_solve(psolver, hm%vhartree, st%rho(1:gr%mesh%np, is))
        folder = 'local.general/potential/'
        write(out_name, '(a,i0,a1,i7.7)')'global-vh.',is,'.',iter
        call dio_function_output(how, &
          trim(folder), trim(out_name), gr%mesh, hm%vhartree, units_out%length, ierr, geo = geo)
      !Computes global XC potential
        call v_ks_calc(ks, hm, st, geo, calc_eigenval = .false. , calc_berry = .false. , calc_energy = .false.)
        folder = 'local.general/potential/'
        write(out_name, '(a,i0,a1,i7.7)')'global-vxc.',is,'.',iter
        call dio_function_output(how, &
          trim(folder), trim(out_name), gr%mesh, hm%vxc(:,is), units_out%length, ierr, geo = geo)
      end do
    end if

    SAFE_DEALLOCATE_A(tmp_vh)
    SAFE_DEALLOCATE_A(st_rho)
    SAFE_DEALLOCATE_A(tmp_rho)
    POP_SUB(local_write_density)
  end subroutine local_write_density

  ! ---------------------------------------------------------
  subroutine local_write_energy(out_energy, nd, domain, lab, inside, center, & 
                                gr, geo, st, hm, ks, mc, iter, l_start, start, how)
    type(local_write_prop_t),      intent(inout) :: out_energy(:)
    integer,                  intent(in)    :: nd 
    type(box_union_t),        intent(in)    :: domain(:)
    character(len=15),        intent(in)    :: lab(:)
    logical,                  intent(in)    :: inside(:,:)
    FLOAT,                    intent(in)    :: center(:,:)
    type(grid_t),         intent(inout) :: gr
    type(geometry_t),     intent(inout) :: geo
    type(states_t),       intent(inout) :: st
    type(hamiltonian_t),  intent(inout) :: hm
    type(v_ks_t),         intent(inout) :: ks
    type(multicomm_t),    intent(in)    :: mc
    integer,              intent(in) :: iter
    integer,              intent(in) :: l_start
    logical,              intent(in) :: start
    integer,              intent(in) :: how

    integer :: id, ii, is, ix, ierr, jd
    character(len=120) :: folder, out_name, aux
    FLOAT              :: geh, gexc, leh, lexc
    FLOAT, allocatable :: tmp_rhoi(:), st_rho(:)
    FLOAT, allocatable :: hm_vxc(:), hm_vh(:)

    PUSH_SUB(local_write_energy)

    SAFE_ALLOCATE(tmp_rhoi(1:gr%mesh%np))
    SAFE_ALLOCATE(st_rho(1:gr%mesh%np_part))
    SAFE_ALLOCATE(hm_vxc(1:gr%mesh%np))
    SAFE_ALLOCATE(hm_vh(1:gr%mesh%np))

    if(.not.mpi_grp_is_root(mpi_world)) return ! only first node outputs

    if(iter == l_start .and. start) then
      do id = 1, nd   
        call local_write_print_header_init(out_energy(id)%handle)

      ! first line -> column names
        write(aux,'(a,i4.4,2x,a,2x,a,2x,a)')'# Domain ID:',id,'- ',trim(lab(id))
        call write_iter_header(out_energy(id)%handle, trim(aux))
        write(aux,'(a)')trim(lab(id))
        call write_iter_header(out_energy(id)%handle, trim(aux))
        call write_iter_nl(out_energy(id)%handle)
        call write_iter_header_start(out_energy(id)%handle)
        aux = 'Total Hartree'
        call write_iter_header(out_energy(id)%handle, trim(aux))
        aux = 'Total Int[n*vxc]'
        call write_iter_header(out_energy(id)%handle, trim(aux))
        write(aux,'(a,i2.2,a,i2.2,a)')'Int[n(',id,')*vh(',id,')]'
        call write_iter_header(out_energy(id)%handle, trim(aux))
        write(aux,'(a,i2.2,a,i2.2,a)')'Int[n(',id,')*vxc(',id,')]'
        call write_iter_header(out_energy(id)%handle, trim(aux))
        do jd = 1, nd
          if (jd /= id) then
            write(aux,'(a,i2.2,a,i2.2,a)')'Int[n(',id,')*vh(',jd,')]'
            call write_iter_header(out_energy(id)%handle, trim(aux))
            write(aux,'(a,i2.2,a,i2.2,a)')'Int[n(',id,')*vxc(',jd,')]'
            call write_iter_header(out_energy(id)%handle, trim(aux))
          end if
        end do
        call write_iter_nl(out_energy(id)%handle)

      ! units
        call write_iter_string(out_energy(id)%handle, '#[Iter n.]')
        call write_iter_header(out_energy(id)%handle, '[' // trim(units_abbrev(units_out%time)) // ']')
        do ii = 1, 2*nd + 2
          call write_iter_header(out_energy(id)%handle, '[' // trim(units_abbrev(units_out%energy)) // ']')
        end do
        call write_iter_nl(out_energy(id)%handle)
      
        call local_write_print_header_end(out_energy(id)%handle)
      end do
    end if

    ! TODO: Make new files for each nspin value. 
    do is = 1, st%d%nspin
      !Compute Hartree potential
      call dpoisson_solve(psolver, hm%vhartree, st%rho(1:gr%mesh%np, is))
      !Compute XC potential
      call v_ks_calc(ks, hm, st, geo, calc_eigenval = .false. , calc_berry = .false. , calc_energy = .false.)
 ! 
      st_rho(:) = st%rho(:, is)
      hm_vxc(:) = hm%vxc(:, is)
      hm_vh(:) = hm%vhartree(:)
     !Compute Global Values
      geh = dmf_integrate(gr%mesh, st_rho(1:gr%mesh%np)*hm%vhartree) 
      gexc = dmf_integrate(gr%mesh, st_rho(1:gr%mesh%np)*hm_vxc(1:gr%mesh%np))
      do id = 1, nd
        call write_iter_set(out_energy(id)%handle, iter)
        call write_iter_start(out_energy(id)%handle)
        call write_iter_double(out_energy(id)%handle, units_from_atomic(units_out%energy, geh), 1)
        call write_iter_double(out_energy(id)%handle, units_from_atomic(units_out%energy, gexc), 1)
     !Compute Local Values
        tmp_rhoi = M_ZERO
        st%rho(:,is) = M_ZERO
        hm%vxc(:,is) = M_ZERO
        hm%vhartree(:) = M_ZERO
        do ix = 1, gr%mesh%np 
          if (inside(ix, id)) st%rho(ix, is) = st_rho(ix)
        end do
        call v_ks_calc(ks, hm, st, geo, calc_eigenval = .false. , calc_berry = .false. , calc_energy = .false.)
        tmp_rhoi(1:gr%mesh%np) = st%rho(1:gr%mesh%np, is)
      !eh = Int[n(id)*v_h(id)]
        leh = dmf_integrate(gr%mesh, tmp_rhoi*hm%vhartree) 
        call write_iter_double(out_energy(id)%handle, units_from_atomic(units_out%energy, leh), 1)
      !exc = Int[n(id)*v_xc(id)]
        lexc = dmf_integrate(gr%mesh, tmp_rhoi*hm%vxc(1:gr%mesh%np,is))
        call write_iter_double(out_energy(id)%handle, units_from_atomic(units_out%energy, lexc), 1)
        do jd = 1, nd
          if (jd /= id) then
            ! TODO: Check for domains overlaping to avoid segmentation fault. 
          !eh = Int[n(id)*v_h(jd)]
            st%rho(:,is) = M_ZERO
            hm%vxc(:,is) = M_ZERO
            hm%vhartree(:) = M_ZERO
            do ix = 1, gr%mesh%np 
              if (inside(ix, jd)) st%rho(ix, is) = st_rho(ix)
            end do
            call v_ks_calc(ks, hm, st, geo, calc_eigenval = .false. , calc_berry = .false. , calc_energy = .false.)
            leh = dmf_integrate(gr%mesh, tmp_rhoi*hm%vhartree) 
            call write_iter_double(out_energy(id)%handle, units_from_atomic(units_out%energy, leh), 1)
          !exc = Int[n(id)*v_xc(jd)]
            lexc = dmf_integrate(gr%mesh, tmp_rhoi*hm%vxc(1:gr%mesh%np, is))
            call write_iter_double(out_energy(id)%handle, units_from_atomic(units_out%energy, lexc), 1)
          end if
        end do
        st%rho(:,is) = st_rho(:)
        hm%vxc(:,is) = hm_vxc(:)
        hm%vhartree(:) = hm_vh(:)
        call write_iter_nl(out_energy(id)%handle)
      end do 
    end do

    SAFE_DEALLOCATE_A(hm_vh)
    SAFE_DEALLOCATE_A(hm_vxc)
    SAFE_DEALLOCATE_A(st_rho)
    SAFE_DEALLOCATE_A(tmp_rhoi)
    POP_SUB(local_write_energy)
  end subroutine local_write_energy

  subroutine local_write_multipole(out_multip, nd, domain, lab, inside, center, & 
                                gr, geo, st, lmax, kick, iter, l_start, start, how)
    type(local_write_prop_t),      intent(inout) :: out_multip(:)
    integer,                  intent(in)    :: nd 
    type(box_union_t),        intent(in)    :: domain(:)
    character(len=15),        intent(in)    :: lab(:)
    logical,                  intent(in)    :: inside(:,:)
    FLOAT,                    intent(in)    :: center(:,:)
    type(grid_t),         intent(in) :: gr
    type(geometry_t),     intent(in) :: geo
    type(states_t),       intent(in) :: st
    integer,              intent(in) :: lmax
    type(kick_t),         intent(in) :: kick
    integer,              intent(in) :: iter
    integer,              intent(in) :: l_start
    logical,              intent(in) :: start
    integer,              intent(in) :: how

    integer :: id, is, ll, mm, add_lm
    character(len=120) :: aux
    FLOAT, allocatable :: ionic_dipole(:,:), multipole(:,:,:)
    CMPLX, allocatable :: zmultipole(:,:,:)
    logical :: cmplxscl

    PUSH_SUB(local_write_multipole)

    cmplxscl = .false.
    if(st%cmplxscl%space) cmplxscl = .true.

    if(mpi_grp_is_root(mpi_world).and. iter == l_start .and. start) then
      do id = 1, nd   
        call local_write_print_header_init(out_multip(id)%handle)
  
        write(aux, '(a15,i2)')      '# nspin        ', st%d%nspin
        call write_iter_string(out_multip(id)%handle, aux)
        call write_iter_nl(out_multip(id)%handle)
  
        write(aux, '(a15,i2)')      '# lmax         ', lmax
        call write_iter_string(out_multip(id)%handle, aux)
        call write_iter_nl(out_multip(id)%handle)

        call kick_write(kick, out = out_multip(id)%handle)

        call write_iter_header_start(out_multip(id)%handle)

        do is = 1, st%d%nspin
          write(aux,'(a18,i1,a1)') 'Electronic charge(', is,')'; call write_iter_header(out_multip(id)%handle, aux)
          if(lmax>0) then
            write(aux, '(a3,a1,i1,a1)') '<x>', '(', is,')'; call write_iter_header(out_multip(id)%handle, aux)
            if(cmplxscl) call write_iter_header(out_multip(id)%handle, ' ')   
            write(aux, '(a3,a1,i1,a1)') '<y>', '(', is,')'; call write_iter_header(out_multip(id)%handle, aux)
            if(cmplxscl) call write_iter_header(out_multip(id)%handle, ' ')   
            write(aux, '(a3,a1,i1,a1)') '<z>', '(', is,')'; call write_iter_header(out_multip(id)%handle, aux)
            if(cmplxscl) call write_iter_header(out_multip(id)%handle, ' ')   
          end if
          do ll = 2, lmax
            do mm = -ll, ll
              write(aux, '(a2,i2,a4,i2,a2,i1,a1)') 'l=', ll, ', m=', mm, ' (', is,')'
              call write_iter_header(out_multip(id)%handle, aux)
            end do
          end do
        end do
        call write_iter_nl(out_multip(id)%handle)

        ! units
        call write_iter_string(out_multip(id)%handle, '#[Iter n.]')
        call write_iter_header(out_multip(id)%handle, '[' // trim(units_abbrev(units_out%time)) // ']')

        do is = 1, st%d%nspin
          do ll = 0, lmax
            do mm = -ll, ll
              select case(ll)
              case(0)
                call write_iter_header(out_multip(id)%handle, 'Electrons')
              case(1)
                call write_iter_header(out_multip(id)%handle, '[' // trim(units_abbrev(units_out%length)) // ']')
                if(cmplxscl) call write_iter_header(out_multip(id)%handle, ' ')   
              case default
                write(aux, '(a,a2,i1)') trim(units_abbrev(units_out%length)), "**", ll
                call write_iter_header(out_multip(id)%handle, '[' // trim(aux) // ']')
                if(cmplxscl) call write_iter_header(out_multip(id)%handle, ' ')   
              end select
            end do
          end do
        end do
        call write_iter_nl(out_multip(id)%handle)

        ! complex quantities
        if(cmplxscl) then
          call write_iter_string(out_multip(id)%handle, '#       _         ')
          call write_iter_header(out_multip(id)%handle, ' ')

          do is = 1, st%d%nspin
            do ll = 0, lmax
              do mm = -ll, ll
                select case(ll)
                case(0)
                  call write_iter_header(out_multip(id)%handle, ' ')
                case(1)
                  call write_iter_header(out_multip(id)%handle, 'Re')
                  call write_iter_header(out_multip(id)%handle, 'Im')   
                case default
                  call write_iter_header(out_multip(id)%handle, 'Re')
                  call write_iter_header(out_multip(id)%handle, 'Im')   
                end select
              end do
            end do
          end do
          call write_iter_nl(out_multip(id)%handle)

        end if
      
        call local_write_print_header_end(out_multip(id)%handle)
        call write_iter_flush(out_multip(id)%handle)
      end do
    end if

    SAFE_ALLOCATE(ionic_dipole(1:gr%mesh%sb%dim, nd))
    SAFE_ALLOCATE(multipole(1:(lmax + 1)**2, 1:st%d%nspin, nd))
    ionic_dipole(:,:) = M_ZERO
    multipole   (:,:,:) = M_ZERO
    if(cmplxscl) then
      SAFE_ALLOCATE(zmultipole(1:(lmax + 1)**2, 1:st%d%nspin, nd))
      zmultipole(:,:,:) = M_z0
    end if

    do is = 1, st%d%nspin
      if(.not. cmplxscl) then
        call dmf_local_multipoles(gr%mesh, nd, st%rho(:,is), lmax, multipole(:,is,:), inside)
      else
        message(1) = 'Local Multipoles is still not implemented for complex densities'
        call messages_fatal(1)
        !FIXME: modify X(mf_local_multipoles) to deal with complex rho
        !call zmf_local_multipoles(gr%mesh, st%zrho%Re(:,is) + M_zI * st%zrho%Im(:,is), lmax,&
        !  zmultipole(:,is), inside, cmplxscl_th = st%cmplxscl%theta, inside)
        !multipole (:,is) = real(zmultipole(:,is)) ! it should be real anyways 
      end if 
    end do
    ! FIXME: with cmplxscl we need to think how to treat 
    ! the ions dipole moment 
    call local_geometry_dipole(nd, domain, geo, ionic_dipole)
    do is = 1, st%d%nspin
      do id = 1, nd
        multipole(2:gr%mesh%sb%dim+1, is, id) = -ionic_dipole(1:gr%mesh%sb%dim, id)/st%d%nspin & 
                                                - multipole(2:gr%mesh%sb%dim+1, is, id)
      end do
    end do

    if(mpi_grp_is_root(mpi_world)) then
      do id = 1, nd
        call write_iter_set(out_multip(id)%handle, iter)
        call write_iter_start(out_multip(id)%handle)
        do is = 1, st%d%nspin
          add_lm = 1
          do ll = 0, lmax
            do mm = -ll, ll
              if(cmplxscl .and. ll > 0 ) then
                message(1) = 'Local Multipoles is still not implemented for complex densities'
                call messages_fatal(1)
                !FIXME: to deal with complex rho
                !call write_iter_double(out_multip(id)%handle, units_from_atomic(units_out%length**ll,&
                !  real(zmultipole(add_lm, is), REAL_PRECISION)), 1)
                !call write_iter_double(out_multip(id)%handle, units_from_atomic(units_out%length**ll,&
                !  aimag(zmultipole(add_lm, is))), 1)
              else
                call write_iter_double(out_multip(id)%handle, units_from_atomic(units_out%length**ll, &
                                        multipole(add_lm, is, id)), 1)
              end if
            add_lm = add_lm + 1
            end do
          end do
        end do
        call write_iter_nl(out_multip(id)%handle)
      end do
    end if

   ! Write multipoles in BILD format
    if(iand(how, OPTION__OUTPUTFORMAT__BILD) /= 0 )then
      !FIXME: to include spin larger than 1.
      is = 1
      do id = 1, nd
        call out_bld_multipoles(multipole(2:4, is, id), center(:,id), lab(id), iter)
      end do
    end if

    SAFE_DEALLOCATE_A(ionic_dipole)
    SAFE_DEALLOCATE_A(multipole)
    SAFE_DEALLOCATE_A(zmultipole)
    POP_SUB(local_write_multipole)
  end subroutine local_write_multipole

  ! ---------------------------------------------------------
  subroutine local_geometry_dipole(nd, dom, geo, dipole)
    integer,           intent(in)  :: nd 
    type(box_union_t), intent(in)  :: dom(:)
    type(geometry_t),  intent(in)  :: geo
    FLOAT,             intent(inout) :: dipole(:,:)

    integer :: ia, id

    PUSH_SUB(local_geometry_dipole)

    dipole(:,:) = M_ZERO
    do ia = 1, geo%natoms
      do  id = 1, nd
        if (box_union_inside(dom(id), geo%atom(ia)%x)) then
          dipole(1:geo%space%dim, id) = dipole(1:geo%space%dim, id) + &
          species_zval(geo%atom(ia)%species)*(geo%atom(ia)%x(1:geo%space%dim))
        end if
      end do
    end do
    dipole = P_PROTON_CHARGE*dipole

    POP_SUB(local_geometry_dipole)
  end subroutine local_geometry_dipole

  ! ---------------------------------------------------------
  subroutine out_bld_multipoles(multipoles, center, label, iter)
    FLOAT,         intent(in) :: multipoles(:)
    FLOAT,         intent(in) :: center(:)
    character(15), intent(in) :: label
    integer,       intent(in) :: iter
   
    integer             :: ll, out_bld
    character(len=80)   :: filename, folder
    FLOAT               :: dipolearrow(3,2)

    PUSH_SUB(out_bld_multipoles)
    
    write(folder,'(a,a)')'local.general/multipoles/',trim(label)
    call io_mkdir(folder)
    write(filename,'(a,a,a,a,i7.7,a)')trim(folder),'/',trim(label),'.',iter,'.bld'
    out_bld = io_open(file=trim(filename), action='write')

    write(out_bld,'(a,a,a,i7)')'.comment ** Arrow for the dipole moment centered at the center of mass for ', &
                        trim(label), ' domain and iteration number: ',iter
    write(out_bld,'(a)')''
    write(out_bld,'(a)')'.color red'
    write(out_bld,'(a,3(f12.6,2x),a)')'.sphere ',(units_from_atomic(units_out%length,center(ll)), ll= 1, 3),' 0.2' 
    do ll = 1, 3
      dipolearrow(ll,1) = units_from_atomic(units_out%length, center(ll) - multipoles(ll))
      dipolearrow(ll,2) = units_from_atomic(units_out%length, center(ll) + multipoles(ll))
    end do
    write(out_bld,'(a,6(f12.6,2x),a)')'.arrow ',(dipolearrow(ll,1), ll= 1, 3), &
                                     (dipolearrow(ll,2), ll= 1, 3), ' 0.1 0.5 0.90'
    call io_close(out_bld)

    POP_SUB(out_bld_multipoles)
  end subroutine out_bld_multipoles

  ! ---------------------------------------------------------
  subroutine local_write_print_header_init(out)
    type(c_ptr), intent(inout) :: out

    PUSH_SUB(local_write_print_header_init)
    call write_iter_clear(out)
    call write_iter_string(out,'################################################################################')
    call write_iter_nl(out)
    call write_iter_string(out,'# HEADER')
    call write_iter_nl(out)

    POP_SUB(local_write_print_header_init)
  end subroutine local_write_print_header_init


  ! ---------------------------------------------------------
  subroutine local_write_print_header_end(out)
    type(c_ptr), intent(inout) :: out

    PUSH_SUB(local_write_print_header_end)

    call write_iter_string(out,'################################################################################')
    call write_iter_nl(out)

    POP_SUB(local_write_print_header_end)
  end subroutine local_write_print_header_end

end module local_write_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
