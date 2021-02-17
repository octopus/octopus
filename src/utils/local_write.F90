!! Copyright (C) 2014 M. Marques, A. Castro, A. Rubio, G. Bertsch, J.Jornet-Somoza
!! Copyright (C) 2020 M. Oliveira
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

module local_write_oct_m
  use iso_c_binding
  use global_oct_m
  use hamiltonian_elec_oct_m
  use io_oct_m
  use io_function_oct_m
  use ions_oct_m
  use kick_oct_m
  use mesh_function_oct_m
  use mesh_oct_m
  use messages_oct_m
  use mpi_oct_m
  use namespace_oct_m
  use parser_oct_m
  use poisson_oct_m
  use profiling_oct_m
  use space_oct_m
  use species_oct_m
  use states_elec_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use varinfo_oct_m
  use v_ks_oct_m
  use write_iter_oct_m

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
    type(local_write_prop_t), allocatable :: out(:)
    integer(8)               :: how              !< how to output
    integer                  :: lmax             !< maximum multipole moment to output
  end type local_write_t

contains

  ! ---------------------------------------------------------
  logical function local_write_check_hm(writ) result(check)
    type(local_write_t), intent(in) :: writ

    PUSH_SUB(local_write_check_hm)

    check = writ%out(LOCAL_OUT_POTENTIAL)%write .or. writ%out(LOCAL_OUT_ENERGY)%write

    PUSH_SUB(local_write_check_hm)
  end function local_write_check_hm

  ! ---------------------------------------------------------
  subroutine local_write_init(writ, namespace, lab, iter, dt)
    type(local_write_t), intent(inout) :: writ
    type(namespace_t),   intent(in)    :: namespace
    character(len=15),   intent(in)    :: lab
    integer,             intent(in)    :: iter
    FLOAT,               intent(in)    :: dt

    integer :: first, flags, iout, default

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
    !%Option energy 8
    !% If set, <tt>octopus</tt> outputs the different components of the energy of the local domains
    !% to the folder <tt>ld.general/energy</tt>.
    !%End

    default = 2**(LOCAL_OUT_MULTIPOLES - 1) 

    call parse_variable(namespace, 'LDOutput', default, flags)

    if (.not.varinfo_valid_option('LDOutput', flags, is_flag = .true.)) then
      call messages_input_error(namespace, 'LDOutput')
    end if

    SAFE_ALLOCATE(writ%out(1:LOCAL_OUT_MAX))
    do iout = 1, LOCAL_OUT_MAX
      writ%out(iout)%write = (bitand(flags, 2**(iout - 1)) /= 0)
    end do

    call messages_obsolete_variable(namespace, 'LDOutputHow', 'LDOutputFormat')
    
    !%Variable LDOutputFormat
    !%Type flag
    !%Default none
    !%Section Utilities::oct-local_multipoles
    !%Description
    !% Describes the format of the output files (see <tt>LDOutput</tt>).
    !% It can take the same values as <tt>OutputFormat</tt> flag.
    !%End
    call parse_variable(namespace, 'LDOutputFormat', 0, writ%how)
    if (.not.varinfo_valid_option('OutputFormat', writ%how, is_flag=.true.)) then
      call messages_input_error(namespace, 'LDOutputFormat')
    end if

    !%Variable LDMultipoleLmax
    !%Type integer
    !%Default 1
    !%Section Utilities::oct-local_multipoles
    !%Description
    !% Maximum electric multipole of the density output to the file <tt>local.multipoles/<>domain%<>.multipoles</tt>
    !% during a time-dependent simulation. Must be non-negative.
    !%End
    call parse_variable(namespace, 'LDMultipoleLmax', 1, writ%lmax)
    if (writ%lmax < 0) then
      write(message(1), '(a,i6,a)') "Input: '", writ%lmax, "' is not a valid LDMultipoleLmax."
      message(2) = '(Must be LDMultipoleLmax >= 0)'
      call messages_fatal(2)
    end if

    call io_mkdir('local.general', namespace)

    if (mpi_grp_is_root(mpi_world)) then
      if (writ%out(LOCAL_OUT_MULTIPOLES)%write) then 
        call io_mkdir('local.general/multipoles', namespace)
        call write_iter_init(writ%out(LOCAL_OUT_MULTIPOLES)%handle, iter, units_from_atomic(units_out%time, dt), &
          trim(io_workpath("local.general/multipoles/"//trim(lab)//".multipoles", namespace)))
      end if

      if (writ%out(LOCAL_OUT_POTENTIAL)%write) then
        call io_mkdir('local.general/potential', namespace)
        call write_iter_init(writ%out(LOCAL_OUT_POTENTIAL)%handle, first, units_from_atomic(units_out%time, dt), &
          trim(io_workpath("local.general/potential/"//trim(lab)//".potential", namespace)))
      end if

      if (writ%out(LOCAL_OUT_DENSITY)%write) then
        call io_mkdir('local.general/densities', namespace)
        call write_iter_init(writ%out(LOCAL_OUT_DENSITY)%handle, first, units_from_atomic(units_out%time, dt), &
          trim(io_workpath("local.general/densities/"//trim(lab)//".densities", namespace)))
      end if

      if (writ%out(LOCAL_OUT_ENERGY)%write) then 
        call io_mkdir('local.general/energy', namespace)
        call write_iter_init(writ%out(LOCAL_OUT_ENERGY)%handle, iter, units_from_atomic(units_out%time, dt), &
          trim(io_workpath("local.general/energy/"//trim(lab)//".energy", namespace)))
      end if
    end if

    POP_SUB(local_write_init)
  end subroutine local_write_init

  ! ---------------------------------------------------------
  subroutine local_write_end(writ)
    type(local_write_t), intent(inout) :: writ

    integer :: iout
    
    PUSH_SUB(local_write_end)

    if (mpi_grp_is_root(mpi_world)) then
      do iout = 1, LOCAL_OUT_MAX
        if (writ%out(iout)%write) then
          call write_iter_end(writ%out(iout)%handle)
        end if
      end do
    end if
    SAFE_DEALLOCATE_A(writ%out)

    POP_SUB(local_write_end)
  end subroutine local_write_end

  ! ---------------------------------------------------------
  subroutine local_write_iter(writ, namespace, space, lab, ions_mask, mesh_mask, mesh, st, hm, ks, ions, kick, iter, l_start, &
                              ldoverwrite)
    type(local_write_t),      intent(inout) :: writ
    type(namespace_t),        intent(in)    :: namespace
    type(space_t),            intent(in)    :: space
    character(len=15),        intent(in)    :: lab
    logical,                  intent(in)    :: ions_mask(:)
    logical,                  intent(in)    :: mesh_mask(:)
    type(mesh_t),             intent(in)    :: mesh
    type(states_elec_t),      intent(inout) :: st
    type(hamiltonian_elec_t), intent(inout) :: hm
    type(v_ks_t),             intent(inout) :: ks
    type(ions_t),             intent(inout) :: ions
    type(kick_t),             intent(inout) :: kick
    integer,                  intent(in)    :: iter
    integer,                  intent(in)    :: l_start
    logical,                  intent(in)    :: ldoverwrite

    type(profile_t), save :: prof

    PUSH_SUB(local_write_iter)
    call profiling_in(prof, "LOCAL_WRITE_ITER")

    if (writ%out(LOCAL_OUT_MULTIPOLES)%write) then
      call local_write_multipole(writ%out(LOCAL_OUT_MULTIPOLES), namespace, lab, ions_mask, mesh_mask, mesh, ions, st, &
        writ%lmax, kick, iter, l_start, ldoverwrite, writ%how)
      if(mpi_grp_is_root(mpi_world)) then
        call write_iter_flush(writ%out(LOCAL_OUT_MULTIPOLES)%handle)
      end if
    end if

    if (writ%out(LOCAL_OUT_DENSITY)%write .or. writ%out(LOCAL_OUT_POTENTIAL)%write) then
      call local_write_density(writ%out(LOCAL_OUT_DENSITY), namespace, space, writ%out(LOCAL_OUT_POTENTIAL), lab, mesh_mask, &
        mesh, ions, st, hm, ks, iter, writ%how)
    end if

    if (writ%out(LOCAL_OUT_ENERGY)%write) then
      call local_write_energy(writ%out(LOCAL_OUT_ENERGY), namespace, space, lab, mesh_mask, mesh, ions, st, hm, ks, iter, l_start, &
        ldoverwrite)
      if(mpi_grp_is_root(mpi_world)) then
        call write_iter_flush(writ%out(LOCAL_OUT_ENERGY)%handle)
      end if
    end if

    call profiling_out(prof)
    POP_SUB(local_write_iter)
  end subroutine local_write_iter

  ! ---------------------------------------------------------
  subroutine local_write_density(out_dens, namespace, space, out_pot, lab, mesh_mask, mesh, ions, st, hm, ks, iter, how)
    type(local_write_prop_t), intent(inout) :: out_dens
    type(namespace_t),        intent(in)    :: namespace
    type(space_t),            intent(in)    :: space
    type(local_write_prop_t), intent(inout) :: out_pot
    character(len=15),        intent(in)    :: lab
    logical,                  intent(in)    :: mesh_mask(:)
    type(mesh_t),             intent(in)    :: mesh
    type(ions_t),             intent(inout) :: ions
    type(states_elec_t),      intent(inout) :: st
    type(hamiltonian_elec_t), intent(inout) :: hm
    type(v_ks_t),             intent(inout) :: ks
    integer,                  intent(in)    :: iter
    integer(8),               intent(in)    :: how

    integer            :: is, ierr
    character(len=120) :: folder, out_name
    FLOAT, allocatable :: tmp_rho(:), st_rho(:)
    FLOAT, allocatable :: tmp_vh(:)

    PUSH_SUB(local_write_density)

    SAFE_ALLOCATE(tmp_rho(1:mesh%np))
    SAFE_ALLOCATE(st_rho(1:mesh%np_part))
    SAFE_ALLOCATE(tmp_vh(1:mesh%np))

    ! FIXME: use just v_ks_calc subroutine for either compute vhartree and vxc
    do is = 1, st%d%nspin
      st_rho(:) = st%rho(:, is)
      where (mesh_mask)
        tmp_rho = st_rho
      elsewhere
        tmp_rho = M_ZERO
      end where

      if (out_dens%write) then
        folder = 'local.general/densities/'//trim(lab)//'.densities/'
        write(out_name, '(a,a1,i0,a1,i7.7)') trim(lab), '.', is,'.', iter
        call dio_function_output(how,  trim(folder), trim(out_name), namespace, space, mesh, tmp_rho(1:mesh%np), &
          units_out%length, ierr, ions = ions)
      end if

      if (out_pot%write) then
        !Computes Hartree potential just for n[r], r belongs to id domain.
        tmp_vh = M_ZERO
        call dpoisson_solve(hm%psolver, tmp_vh, tmp_rho)
        folder = 'local.general/potential/'//trim(lab)//'.potential/'
        write(out_name, '(a,i0,a1,i7.7)') 'vh.', is, '.', iter
        call dio_function_output(how, trim(folder), trim(out_name), namespace, space, mesh, tmp_vh, units_out%length, ierr, &
          ions = ions)
        !Computes XC potential
        where (mesh_mask)
          st%rho(:, is) = st_rho
        elsewhere
          st%rho(:, is) = M_ZERO
        end where
        call v_ks_calc(ks, namespace, space, hm, st, ions, calc_eigenval = .false. , calc_energy = .false.)
        folder = 'local.general/potential/'//trim(lab)//'.potential/'
        write(out_name, '(a,i0,a1,i7.7)') 'vxc.', is, '.', iter
        call dio_function_output(how, trim(folder), trim(out_name), namespace, space, mesh, hm%vxc(:,is), units_out%length, ierr, &
          ions = ions)
        st%rho(:,is) = st_rho(:)
      end if
    end do

    if (out_pot%write) then
      do is = 1, st%d%nspin
        !Computes Hartree potential
        call dpoisson_solve(hm%psolver, hm%vhartree, st%rho(1:mesh%np, is))
        folder = 'local.general/potential/'
        write(out_name, '(a,i0,a1,i7.7)') 'global-vh.', is, '.', iter
        call dio_function_output(how, trim(folder), trim(out_name), namespace, space, mesh, hm%vhartree, units_out%length, ierr, &
          ions = ions)
        !Computes global XC potential
        call v_ks_calc(ks, namespace, space, hm, st, ions, calc_eigenval = .false. , calc_energy = .false.)
        folder = 'local.general/potential/'
        write(out_name, '(a,i0,a1,i7.7)') 'global-vxc.', is, '.', iter
        call dio_function_output(how, trim(folder), trim(out_name), namespace, space, mesh, hm%vxc(:,is), units_out%length, ierr, &
          ions = ions)
      end do
    end if

    SAFE_DEALLOCATE_A(tmp_vh)
    SAFE_DEALLOCATE_A(st_rho)
    SAFE_DEALLOCATE_A(tmp_rho)

    POP_SUB(local_write_density)
  end subroutine local_write_density

  ! ---------------------------------------------------------
  subroutine local_write_energy(out_energy, namespace, space, lab, mesh_mask, mesh, ions, st, hm, ks, iter, l_start, start)
    type(local_write_prop_t), intent(inout) :: out_energy
    type(namespace_t),        intent(in)    :: namespace
    type(space_t),            intent(in)    :: space
    character(len=15),        intent(in)    :: lab
    logical,                  intent(in)    :: mesh_mask(:)
    type(mesh_t),             intent(in)    :: mesh
    type(ions_t),             intent(inout) :: ions
    type(states_elec_t),      intent(inout) :: st
    type(hamiltonian_elec_t), intent(inout) :: hm
    type(v_ks_t),             intent(inout) :: ks
    integer,                  intent(in)    :: iter
    integer,                  intent(in)    :: l_start
    logical,                  intent(in)    :: start

    !integer :: ix
    integer :: ii, is
    character(len=120) :: aux
    FLOAT              :: geh, gexc, leh, lexc
    FLOAT, allocatable :: tmp_rhoi(:), st_rho(:)
    FLOAT, allocatable :: hm_vxc(:), hm_vh(:)

    PUSH_SUB(local_write_energy)

    SAFE_ALLOCATE(tmp_rhoi(1:mesh%np))
    SAFE_ALLOCATE(st_rho(1:mesh%np_part))
    SAFE_ALLOCATE(hm_vxc(1:mesh%np))
    SAFE_ALLOCATE(hm_vh(1:mesh%np_part))

    if(mpi_grp_is_root(mpi_world)) then  ! only first node outputs

      if(iter == l_start .and. start) then
        call local_write_print_header_init(out_energy%handle)
    
        ! first line -> column names
        write(aux,'(a,2x,a)')'# Domain ID:', trim(lab)
        call write_iter_header(out_energy%handle, trim(aux))
        write(aux,'(a)') trim(lab)
        call write_iter_header(out_energy%handle, trim(aux))
        call write_iter_nl(out_energy%handle)
        call write_iter_header_start(out_energy%handle)
        aux = 'Total Hartree'
        call write_iter_header(out_energy%handle, trim(aux))
        aux = 'Total Int[n*vxc]'
        call write_iter_header(out_energy%handle, trim(aux))
        write(aux,'(a)') 'Int[n*vh]'
        call write_iter_header(out_energy%handle, trim(aux))
        write(aux,'(a)')'Int[n*vxc]'
        call write_iter_header(out_energy%handle, trim(aux))
! Cross-terms between different domains are currently disabled
!        do jd = 1, nd
!          if (jd /= id) then
!            write(aux,'(a,i2.2,a,i2.2,a)')'Int[n(',id,')*vh(',jd,')]'
!            call write_iter_header(out_energy(id)%handle, trim(aux))
!            write(aux,'(a,i2.2,a,i2.2,a)')'Int[n(',id,')*vxc(',jd,')]'
!            call write_iter_header(out_energy(id)%handle, trim(aux))
!          end if
!        end do
        call write_iter_nl(out_energy%handle)
    
        ! units
        call write_iter_string(out_energy%handle, '#[Iter n.]')
        call write_iter_header(out_energy%handle, '[' // trim(units_abbrev(units_out%time)) // ']')
        do ii = 1, 4
          call write_iter_header(out_energy%handle, '[' // trim(units_abbrev(units_out%energy)) // ']')
        end do
        call write_iter_nl(out_energy%handle)
        
        call local_write_print_header_end(out_energy%handle)
      end if
    end if

    ! TODO: Make new files for each nspin value. 
    do is = 1, st%d%nspin
      ! Compute Hartree potential
      call dpoisson_solve(hm%psolver, hm%vhartree, st%rho(1:mesh%np, is))
      ! Compute XC potential
      call v_ks_calc(ks, namespace, space, hm, st, ions, calc_eigenval = .false. , calc_energy = .false.)

      st_rho(:) = st%rho(:, is)
      hm_vxc(:) = hm%vxc(:, is)
      hm_vh(:) = hm%vhartree(:)
      ! Compute Global Values
      geh = dmf_integrate(mesh, st_rho(1:mesh%np)*hm%vhartree(1:mesh%np)) 
      gexc = dmf_integrate(mesh, st_rho(1:mesh%np)*hm_vxc(1:mesh%np))
      if (mpi_grp_is_root(mpi_world)) then
        call write_iter_set(out_energy%handle, iter)
        call write_iter_start(out_energy%handle)
        call write_iter_double(out_energy%handle, units_from_atomic(units_out%energy, geh), 1)
        call write_iter_double(out_energy%handle, units_from_atomic(units_out%energy, gexc), 1)
      end if

      !Compute Local Values
      hm%vxc(:,is) = M_ZERO
      hm%vhartree(:) = M_ZERO
      where (mesh_mask)
        st%rho(:, is) = st_rho
      elsewhere
        st%rho(:, is) = M_ZERO
      end where

      call v_ks_calc(ks, namespace, space, hm, st, ions, calc_eigenval = .false. , calc_energy = .false.)
      tmp_rhoi(1:mesh%np) = st%rho(1:mesh%np, is)
      !eh = Int[n*v_h]
      leh = dmf_integrate(mesh, tmp_rhoi*hm%vhartree(1:mesh%np)) 
      !exc = Int[n*v_xc]
      lexc = dmf_integrate(mesh, tmp_rhoi*hm%vxc(1:mesh%np,is))
      if (mpi_grp_is_root(mpi_world)) then
        call write_iter_double(out_energy%handle, units_from_atomic(units_out%energy, leh), 1)
        call write_iter_double(out_energy%handle, units_from_atomic(units_out%energy, lexc), 1)
      end if

! Cross-terms between different domains are currently disabled
!      do jd = 1, nd
!        if (jd /= id) then
          ! TODO: Check for domains overlaping to avoid segmentation fault. 
          !eh = Int[n(id)*v_h(jd)]
!          st%rho(:,is) = M_ZERO
!          hm%vxc(:,is) = M_ZERO
!          hm%vhartree(:) = M_ZERO
!          do ix = 1, mesh%np 
!            if (mesh_mask(ix, jd)) st%rho(ix, is) = st_rho(ix)
!          end do
!          call v_ks_calc(ks, namespace, hm, st, ions, calc_eigenval = .false. , calc_energy = .false.)
!          leh = dmf_integrate(mesh, tmp_rhoi*hm%vhartree(1:mesh%np)) 
          !exc = Int[n(id)*v_xc(jd)]
!          lexc = dmf_integrate(mesh, tmp_rhoi*hm%vxc(1:mesh%np, is))
!          if (mpi_grp_is_root(mpi_world)) then
!            call write_iter_double(out_energy(id)%handle, units_from_atomic(units_out%energy, leh), 1)
!            call write_iter_double(out_energy(id)%handle, units_from_atomic(units_out%energy, lexc), 1)
!          end if
!        end if
!      end do
      st%rho(:,is) = st_rho(:)
      hm%vxc(:,is) = hm_vxc(:)
      hm%vhartree(:) = hm_vh(:)
      if (mpi_grp_is_root(mpi_world)) then
        call write_iter_nl(out_energy%handle)
      end if
    end do

    SAFE_DEALLOCATE_A(hm_vh)
    SAFE_DEALLOCATE_A(hm_vxc)
    SAFE_DEALLOCATE_A(st_rho)
    SAFE_DEALLOCATE_A(tmp_rhoi)

    POP_SUB(local_write_energy)
  end subroutine local_write_energy

  ! ---------------------------------------------------------
  subroutine local_write_multipole(out_multip, namespace, lab, ions_mask, mesh_mask, mesh, ions, st, lmax, kick, iter, l_start, &
    start, how)
    type(local_write_prop_t), intent(inout) :: out_multip
    type(namespace_t),        intent(in)    :: namespace
    character(len=15),        intent(in)    :: lab
    logical,                  intent(in)    :: ions_mask(:)
    logical,                  intent(in)    :: mesh_mask(:)
    type(mesh_t),             intent(in)    :: mesh
    type(ions_t),             intent(in)    :: ions
    type(states_elec_t),      intent(in)    :: st
    integer,                  intent(in)    :: lmax
    type(kick_t),             intent(in)    :: kick
    integer,                  intent(in)    :: iter
    integer,                  intent(in)    :: l_start
    logical,                  intent(in)    :: start
    integer(8),               intent(in)    :: how

    integer :: is, ll, mm, add_lm
    character(len=120) :: aux
    FLOAT :: ionic_dipole(ions%space%dim), center(MAX_DIM)
    FLOAT, allocatable :: multipole(:,:)
    logical :: use_ionic_dipole

    PUSH_SUB(local_write_multipole)

    if (mpi_grp_is_root(mpi_world).and. iter == l_start .and. start) then
      call local_write_print_header_init(out_multip%handle)
  
      write(aux, '(a15,i2)')      '# nspin        ', st%d%nspin
      call write_iter_string(out_multip%handle, aux)
      call write_iter_nl(out_multip%handle)
  
      write(aux, '(a15,i2)')      '# lmax         ', lmax
      call write_iter_string(out_multip%handle, aux)
      call write_iter_nl(out_multip%handle)

      call kick_write(kick, out = out_multip%handle)

      call write_iter_header_start(out_multip%handle)

      do is = 1, st%d%nspin
        write(aux,'(a18,i1,a1)') 'Electronic charge(', is,')'; call write_iter_header(out_multip%handle, aux)
        if (lmax > 0) then
          write(aux, '(a3,a1,i1,a1)') '<x>', '(', is,')'; call write_iter_header(out_multip%handle, aux)
          write(aux, '(a3,a1,i1,a1)') '<y>', '(', is,')'; call write_iter_header(out_multip%handle, aux)
          write(aux, '(a3,a1,i1,a1)') '<z>', '(', is,')'; call write_iter_header(out_multip%handle, aux)
        end if
        do ll = 2, lmax
          do mm = -ll, ll
            write(aux, '(a2,i2,a4,i2,a2,i1,a1)') 'l=', ll, ', m=', mm, ' (', is,')'
            call write_iter_header(out_multip%handle, aux)
          end do
        end do
      end do
      call write_iter_nl(out_multip%handle)

      ! units
      call write_iter_string(out_multip%handle, '#[Iter n.]')
      call write_iter_header(out_multip%handle, '[' // trim(units_abbrev(units_out%time)) // ']')

      do is = 1, st%d%nspin
        do ll = 0, lmax
          do mm = -ll, ll
            select case(ll)
            case(0)
              call write_iter_header(out_multip%handle, 'Electrons')
            case(1)
              call write_iter_header(out_multip%handle, '[' // trim(units_abbrev(units_out%length)) // ']')
            case default
              write(aux, '(a,a2,i1)') trim(units_abbrev(units_out%length)), "**", ll
              call write_iter_header(out_multip%handle, '[' // trim(aux) // ']')
            end select
          end do
        end do
      end do
      call write_iter_nl(out_multip%handle)

      call local_write_print_header_end(out_multip%handle)
      call write_iter_flush(out_multip%handle)
    end if

    center(1:ions%space%dim) = ions%center_of_mass(mask=ions_mask)

    SAFE_ALLOCATE(multipole(1:(lmax + 1)**2, 1:st%d%nspin))
    multipole = M_ZERO

    do is = 1, st%d%nspin
      call dmf_multipoles(mesh, st%rho(:,is), lmax, multipole(:,is), mask=mesh_mask)
    end do
    ! Setting center of mass as reference needed for non-neutral systems.
    do is = 1, st%d%nspin
      multipole(2:mesh%sb%dim+1, is) =  multipole(2:mesh%sb%dim+1, is) - center(1:mesh%sb%dim)*multipole(1, is)
    end do

    ! For transition densities or density differences there is no need to add the ionic dipole
    
    !%Variable LDIonicDipole
    !%Type logical
    !%Default yes
    !%Section Utilities::oct-local_multipoles
    !%Description
    !% Describes if the ionic dipole has to be take into account
    !% when computing the multipoles.
    !%End
    call parse_variable(namespace, 'LDIonicDipole', .true., use_ionic_dipole)
    if (use_ionic_dipole) then
      ! Calculate ionic dipole, setting the center of mass as reference, as that is needed for non-neutral systems.
      ionic_dipole(1:ions%space%dim) = ions%dipole(mask=ions_mask) + &
        P_PROTON_CHARGE*ions%val_charge(mask=ions_mask)*center(1:ions%space%dim)

      do is = 1, st%d%nspin
        multipole(2:ions%space%dim+1, is) = -ionic_dipole(1:ions%space%dim)/st%d%nspin - multipole(2:ions%space%dim+1, is)
      end do
    end if

    if (mpi_grp_is_root(mpi_world)) then
      call write_iter_set(out_multip%handle, iter)
      call write_iter_start(out_multip%handle)
      do is = 1, st%d%nspin
        add_lm = 1
        do ll = 0, lmax
          do mm = -ll, ll
            call write_iter_double(out_multip%handle, units_from_atomic(units_out%length**ll, multipole(add_lm, is)), 1)
            add_lm = add_lm + 1
          end do
        end do
      end do
      call write_iter_nl(out_multip%handle)
    end if

    ! Write multipoles in BILD format
    if(bitand(how, OPTION__OUTPUTFORMAT__BILD) /= 0 .and. mpi_grp_is_root(mpi_world))then
      !FIXME: to include spin larger than 1.
      is = 1
      call out_bld_multipoles(namespace, multipole(2:4, is), center, lab, iter)
    end if

    SAFE_DEALLOCATE_A(multipole)

    POP_SUB(local_write_multipole)
  end subroutine local_write_multipole

  ! ---------------------------------------------------------
  subroutine out_bld_multipoles(namespace, multipoles, center, label, iter)
    type(namespace_t), intent(in) :: namespace
    FLOAT,             intent(in) :: multipoles(:)
    FLOAT,             intent(in) :: center(:)
    character(15),     intent(in) :: label
    integer,           intent(in) :: iter
   
    integer             :: ll, out_bld
    character(len=80)   :: filename, folder
    FLOAT               :: dipolearrow(3,2)

    PUSH_SUB(out_bld_multipoles)
    
    write(folder,'(a,a)')'local.general/multipoles/',trim(label)
    call io_mkdir(folder, namespace)
    write(filename,'(a,a,a,a,i7.7,a)')trim(folder),'/',trim(label),'.',iter,'.bld'
    out_bld = io_open(trim(filename), namespace, action='write')

    write(out_bld,'(a,a,a,i7)')'.comment ** Arrow for the dipole moment centered at the center of mass for ', &
                        trim(label), ' domain and iteration number: ',iter
    write(out_bld,'(a)')''
    write(out_bld,'(a)')'.color red'
    write(out_bld,'(a,3(f12.6,2x),a)')'.sphere ',(units_from_atomic(units_out%length,center(ll)), ll= 1, 3),' 0.2' 
    do ll = 1, 3
      dipolearrow(ll,1) = units_from_atomic(units_out%length, center(ll) - multipoles(ll)/M_TWO)
      dipolearrow(ll,2) = units_from_atomic(units_out%length, center(ll) + multipoles(ll)/M_TWO)
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

end module local_write_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
