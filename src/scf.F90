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

#include "global.h"

module scf
  use global
  use lib_oct_parser
  use units
  use geometry
  use states
  use hamiltonian
  use eigen_solver
  use mix
  use lcao

  implicit none

  integer, parameter :: MIXDENS = 0, &
                        MIXPOT  = 1

  type scf_type  ! some variables used for the scf cycle
    integer :: max_iter ! maximum number of scf iterations

    FLOAT :: conv_abs_dens, conv_rel_dens, &
             conv_abs_ev, conv_rel_ev         ! several convergence criteria
    FLOAT :: abs_dens, rel_dens, abs_ev, rel_ev

    integer :: what2mix

    logical :: lcao_restricted

    type(mix_type) :: smix
    type(eigen_solver_type) :: eigens
  end type scf_type

contains

subroutine scf_init(scf, m, st, h)
  type(scf_type),         intent(inout) :: scf
  type(mesh_type),        intent(in)    :: m
  type(states_type),      intent(in)    :: st
  type(hamiltonian_type), intent(in)    :: h

  call push_sub('scf_init')

  call loct_parse_int  ("MaximumIter",        200, scf%max_iter)
  call loct_parse_float("ConvAbsDens", CNST(1e-5), scf%conv_abs_dens)
  call loct_parse_float("ConvRelDens",     M_ZERO, scf%conv_rel_dens)
  call loct_parse_float("ConvAbsEv",       M_ZERO, scf%conv_abs_ev)
  call loct_parse_float("ConvRelEv",       M_ZERO, scf%conv_rel_ev)

  if(scf%max_iter <= 0 .and. &
      scf%conv_abs_dens <= M_ZERO .and. scf%conv_rel_dens <= M_ZERO .and. &
      scf%conv_abs_ev <= M_ZERO .and. scf%conv_rel_ev <= M_ZERO) then
    message(1) = "Input: Not all convergence criteria can be <= 0"
    message(2) = "Please set one of the following:"
    message(3) = "MaximumIter | ConvAbsDens | ConvRelDens | ConvAbsEv | ConvRelEv"
    call write_fatal(3)
  end if

  if(scf%max_iter <= 0) scf%max_iter = huge(scf%max_iter)

  call loct_parse_int("What2Mix", 0, scf%what2mix)
  select case (scf%what2mix)
  case (MIXDENS)
     message(1) = 'Info: SCF mixing the density.'
  case (MIXPOT)
     if (h%ip_app) then
        message(1) = "Input: Cannot mix the potential whit non-interacting electrons."
        call write_fatal(1)
     else
        message(1) = 'Info: SCF mixing the potential.'
     end if
  case default
      write(message(1), '(a,i5,a)') "Input: '", scf%what2mix, &
           "' is not a valid option for What2Mix"
      message(2) = '(What2Mix = 0 | 1)'
      call write_fatal(2)
   end select
   call write_info(1)

  ! Handle mixing now...
  call mix_init(scf%smix, m%np, st%d%nspin)

  ! now the eigen solver stuff
  call eigen_solver_init(scf%eigens, st, m)

  ! Should the calculation be restricted to LCAO subspace?
  call loct_parse_logical("SCFinLCAO", .false., scf%lcao_restricted)
  if(scf%lcao_restricted) then
    message(1) = 'Info: SCF restricted to LCAO subspace'
    call write_info(1)
  endif

  call pop_sub(); return
end subroutine scf_init

subroutine scf_end(scf)
  type(scf_type), intent(inout) :: scf

  call push_sub('scf_end') 

  call eigen_solver_end(scf%eigens)
  call mix_end(scf%smix)

  call pop_sub(); return
end subroutine scf_end

subroutine scf_run(scf, m, st, geo, h, outp)
  type(scf_type),         intent(inout) :: scf
  type(mesh_type),        intent(in)    :: m
  type(states_type),      intent(inout) :: st
  type(geometry_type),    intent(inout) :: geo
  type(hamiltonian_type), intent(inout) :: h
  type(output_type),      intent(in)    :: outp

  integer :: iter, iunit, is
  FLOAT :: evsum_out, evsum_in
  FLOAT, allocatable :: rhoout(:,:), rhoin(:,:)
  FLOAT, allocatable :: vout(:,:), vin(:,:)
  logical :: finish

  call push_sub('scf_run')

  if(scf%lcao_restricted) call lcao_init(m, st, geo, h)

  allocate(rhoout(m%np, st%d%nspin), rhoin(m%np, st%d%nspin))
  rhoin = st%rho
  if (scf%what2mix == MIXPOT) then
     allocate(vout(m%np, st%d%nspin), vin(m%np, st%d%nspin))
     vin = h%vhxc; vout = M_ZERO
  end if
  evsum_in = states_eigenvalues_sum(st)

  ! SCF cycle
  do iter = 1, scf%max_iter
    if(scf%lcao_restricted) then
      call lcao_wf(m, st, h)
    else
      call eigen_solver_run(scf%eigens, m, st, h, iter)
    endif

    ! occupations
    call states_fermi(st, m)

    ! compute output density, potential (if needed) and eigenvalues sum
    call X(calcdens)(st, m%np, rhoout)
    if (scf%what2mix == MIXPOT) then
       st%rho = rhoout
       call X(h_calc_vhxc) (h, m, st)
       vout = h%vhxc
    end if
    evsum_out = states_eigenvalues_sum(st)

    ! compute convergence criteria
    scf%abs_dens = M_ZERO
    do is = 1, st%d%nspin
       scf%abs_dens = scf%abs_dens + dmf_integrate(m, (rhoin(:,is) - rhoout(:,is))**2)
    end do
    scf%abs_dens = sqrt(scf%abs_dens)
    scf%rel_dens = scf%abs_dens / st%qtot
    scf%abs_ev = abs(evsum_out - evsum_in)
    scf%rel_ev = scf%abs_ev / abs(evsum_out)

    ! are we finished?
    finish = &
        (scf%conv_abs_dens > M_ZERO .and. scf%abs_dens <= scf%conv_abs_dens) .or. &
        (scf%conv_rel_dens > M_ZERO .and. scf%rel_dens <= scf%conv_rel_dens) .or. &
        (scf%conv_abs_ev   > M_ZERO .and. scf%abs_ev   <= scf%conv_abs_ev) .or. &
        (scf%conv_rel_ev   > M_ZERO .and. scf%rel_ev   <= scf%conv_rel_ev)

    call scf_write_iter

    select case (scf%what2mix)
    case (MIXDENS)
       ! mix input and output densities and compute new potential
       call mixing(scf%smix, iter, m%np, st%d%nspin, rhoin, rhoout, st%rho)
       call X(h_calc_vhxc) (h, m, st)
    case (MIXPOT)
       ! mix input and output potentials
       call mixing(scf%smix, iter, m%np, st%d%nspin, vin, vout, h%vhxc)
    end select

    ! save restart information
    if(finish.or.(modulo(iter, 3) == 0).or.iter==scf%max_iter.or.clean_stop()) &
         call X(states_write_restart)("tmp/restart.static", m, st)

    if(finish) then
      write(message(1), '(a, i4, a)')'Info: SCF converged in ', iter, ' iterations'
      call write_info(1)
      if(scf%lcao_restricted) call lcao_end
      exit
    end if

    if(outp%duringscf) then
      call X(states_output) (st, m, "static", outp)
      call hamiltonian_output(h, m, "static", outp)
    endif

    ! save information for the next iteration
    rhoin = st%rho
    if (scf%what2mix == MIXPOT) vin = h%vhxc
    evsum_in = evsum_out

    if(clean_stop()) exit
  end do

  if (scf%what2mix == MIXPOT) then
     call X(h_calc_vhxc) (h, m, st)
     deallocate(vout, vin)
  end if
  deallocate(rhoout, rhoin)

  if(.not.finish) then
    message(1) = 'SCF *not* converged!'
    call write_warning(1)
  end if

  ! calculate forces
  call X(epot_forces)(h%ep, m, st, geo)

  ! output final information
  call scf_write_static("static", "info")
  call X(states_output) (st, m, "static", outp)
  if(outp%what(output_geometry)) &
     call atom_write_xyz("static", "geometry", geo)
  call hamiltonian_output(h, m, "static", outp)

  if (conf%periodic_dim>0.and.st%nik>st%d%nspin) then
    call io_assign(iunit)
    open(iunit, status='unknown', file='static/bands.dat')
    call states_write_bands(iunit, st%nst, st)
    call io_close(iunit)
  end if

  call pop_sub()
contains

subroutine scf_write_iter
  write(message(1),'(a)') '************'
  write(message(2),'(a,i5)') 'SCF CYCLE ITER #',iter
  write(message(3),'(2(a,es9.2))') &
         ' abs_dens = ', scf%abs_dens, ' abs_ev = ', scf%abs_ev
  write(message(4),'(2(a,es9.2))') &
         ' rel_dens = ', scf%rel_dens, ' rel_ev = ', scf%rel_ev
  call write_info(4)
  if(.not.scf%lcao_restricted) then
    write(message(1),'(a,i6)') 'Matrix vector products: ', scf%eigens%matvec
    write(message(2),'(a,i6)') 'Converged eigenvectors: ', scf%eigens%converged
    call write_info(2)
    call states_write_eigenvalues(stdout, st%nst, st, scf%eigens%diff)
  else
    call states_write_eigenvalues(stdout, st%nst, st)
  endif
  write(message(1),'(a)') '************'
  write(message(2),'(a)')
  call write_info(2)
end subroutine scf_write_iter

subroutine scf_write_static(dir, fname)
  character(len=*), intent(in) :: dir, fname

  FLOAT :: e_dip(conf%dim, st%d%nspin), n_dip(conf%dim)
  FLOAT, parameter :: ATOMIC_TO_DEBYE = CNST(2.5417462)
  integer :: iunit, i, j

  call loct_mkdir(trim(dir))
  call io_assign(iunit)
  open(iunit, status='unknown', file=trim(dir) // "/" // trim(fname))

  ! mesh
  write(iunit, '(a,a)') 'System name: ', geo%sysname
  write(iunit, '(1x)')

  write(iunit, '(a)') 'Mesh:'
  call mesh_write_info(m, iunit)
  write(iunit,'(1x)')
  
  if (conf%periodic_dim > 0) then
    call kpoints_write_info(st,iunit)
    write(iunit,'(1x)')
  end if

  if(.not. h%ip_app) then
    write(iunit, '(a)') 'Exchange and correlation functionals:'
    call xc_write_info(h%xc, iunit)
  else
    write(iunit, '(a)') 'Independent Particles'
  end if
  write(iunit,'(1x)')

  ! scf information
  if(finish) then
    write(iunit, '(a, i4, a)')'SCF converged in ', iter, ' iterations'
  else
    write(iunit, '(a)') 'SCF *not* converged!'
  end if
  write(iunit, '(1x)')

  call states_write_eigenvalues(iunit, st%nst, st)
  write(iunit, '(1x)')

  write(iunit, '(a)') 'Energy:'
  call hamiltonian_energy(h, st, geo%eii, iunit)
  write(iunit, '(1x)')

  if(st%d%ispin > UNPOLARIZED) then
    call write_magnet(iunit, m, st)
  end if

  ! Next lines of code calculate the dipole of the molecule, summing the electronic and
  ! ionic contributions.
                 call states_calculate_multipoles(m, st, (/ M_ONE, M_ZERO, M_ZERO /), e_dip(1, :))
  if(conf%dim>1) call states_calculate_multipoles(m, st, (/ M_ZERO, M_ONE, M_ZERO /), e_dip(2, :))
  if(conf%dim>2) call states_calculate_multipoles(m, st, (/ M_ZERO, M_ZERO, M_ONE /), e_dip(3, :))
  do j = 1, conf%dim
     e_dip(j, 1) = sum(e_dip(j, :))
  enddo
  call geometry_dipole(geo, n_dip)
  n_dip(:) = n_dip(:) - e_dip(:, 1)
  write(iunit, '(3a)') 'Dipole [', trim(units_out%length%abbrev), ']:                    [Debye]'
  do j = 1, conf%dim
     write(iunit, '(6x,a,i1,a,es14.5,3x,2es14.5)') '<x', j, '> = ', n_dip(j) / units_out%length%factor, &
                                                                    n_dip(j)*ATOMIC_TO_DEBYE
  enddo
  write(iunit,'(a)')

  write(iunit, '(a)') 'Convergence:'
  write(iunit, '(6x, a, es14.8,a,es14.8,a)') 'abs_dens = ', scf%abs_dens, &
      ' (', scf%conv_abs_dens, ')'
  write(iunit, '(6x, a, es14.8,a,es14.8,a)') 'rel_dens = ', scf%rel_dens, &
      ' (', scf%conv_rel_dens, ')'
  write(iunit, '(6x, a, es14.8,a,es14.8,4a)') 'abs_ev = ', scf%abs_ev, &
      ' (', scf%conv_abs_ev / units_out%energy%factor, ')', &
      ' [',  trim(units_out%energy%abbrev), ']'
  write(iunit, '(6x, a, es14.8,a,es14.8,a)') 'rel_ev = ', scf%rel_ev, &
      ' (', scf%conv_rel_ev, ')'
  write(iunit,'(1x)') 

  write(iunit,'(3a)') 'Forces on the ions [', trim(units_out%force%abbrev), "]"
  write(iunit,'(a,10x,14x,a,14x,a,14x,a)') ' Ion','x','y','z'
  do i = 1,geo%natoms
    write(iunit,'(i4,a10,3f15.6)') i, trim(geo%atom(i)%spec%label), &
         geo%atom(i)%f(:) / units_out%force%factor
  end do

  call io_close(iunit)
end subroutine scf_write_static

subroutine write_magnet(iunit, mesh, st)
  integer, intent(in) :: iunit
  type(mesh_type), intent(in) :: mesh
  type(states_type), intent(IN) :: st
  
  FLOAT :: m(3), sign
  R_TYPE :: c
  integer :: i, ik, ist
  
  write(iunit, '(a)') 'Magnetization:'
  if(st%d%ispin == SPIN_POLARIZED) then ! collinear spin
    sign = M_ONE
    m(3) = M_ZERO
    do ik = 1, st%nik
      do ist = 1, st%nst
        m(3) = m(3) + sign*st%d%kweights(ik)*st%occ(ist, ik)
      end do
      sign = -sign
    end do
    write(iunit, '(a,f15.6)') ' mz = ', m(3)
    
  else if(st%d%ispin == SPINORS) then ! non-collinear
    m = M_ZERO
    do ik = 1, st%nik
      do ist = 1, st%nst
        do i = 1, mesh%np
          c = R_CONJ(st%X(psi) (i, 1, ist, ik)) * st%X(psi) (i, 2, ist, ik)
          m(1) = m(1) + st%d%kweights(ik)*st%occ(ist, ik)* M_TWO*R_REAL(c)
          m(2) = m(2) + st%d%kweights(ik)*st%occ(ist, ik)* M_TWO*R_AIMAG(c)
          c = R_ABS(st%X(psi) (i, 1, ist, ik))**2 - R_ABS(st%X(psi) (i, 2, ist, ik))**2
          m(3) = m(3) + st%d%kweights(ik)*st%occ(ist, ik)* R_REAL(c)
        end do
      end do
    end do
    m = m*mesh%vol_pp
    write(iunit, '(a,f15.6)') ' mx = ', m(1)
    write(iunit, '(a,f15.6)') ' my = ', m(2)
    write(iunit, '(a,f15.6)') ' mz = ', m(3)
  end if
  
  write(iunit,'(1x)') 
  
end subroutine write_magnet

end subroutine scf_run

end module scf
