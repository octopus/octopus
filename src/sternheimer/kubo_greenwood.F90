!! Copyright (C) 2018 X. Andrade
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

module kubo_greenwood_oct_m
  use boundaries_oct_m
  use derivatives_oct_m
  use global_oct_m
  use hamiltonian_oct_m
  use io_oct_m
  use messages_oct_m
  use mesh_oct_m
  use mesh_function_oct_m
  use restart_oct_m
  use smear_oct_m
  use states_oct_m
  use states_restart_oct_m
  use system_oct_m
  use unit_oct_m
  use unit_system_oct_m
  use parser_oct_m
  use profiling_oct_m
 
  implicit none

  private
  public ::                       &
       kubo_greenwood_run

contains
  
  subroutine kubo_greenwood_run(sys, hm)
    type(system_t),       intent(inout) :: sys
    type(hamiltonian_t),  intent(inout) :: hm

    type(restart_t) :: gs_restart
    integer :: ierr, nfreq, ifreq, iunit
    integer :: ist, jst, iqn, idim, idir, jdir
    CMPLX, allocatable :: psii(:, :), psij(:, :), gpsii(:, :, :), gpsij(:, :, :)
    CMPLX, allocatable :: tensor(:, :, :), therm_tensor(:,:,:)
    CMPLX, allocatable :: k11(:,:,:), k12(:,:,:), k21(:,:,:), k22(:,:,:)
    CMPLX :: prod
    FLOAT :: eigi, eigj, occi, occj, df, width, dfreq, maxfreq, minfreq, ww
    FLOAT :: sum_occ, sum_cond1, sum_cond2, sum_cond, sum_k, e_fermi
    type(mesh_t), pointer :: mesh
    character(len=80) :: dirname, str_tmp

    PUSH_SUB(kubo_greewood_run)

    call messages_write('Info: Starting Kubo-Greenwood linear-response calculation.')
    call messages_info()
    
    call messages_experimental('Kubo Greenwood')

      !%Variable KuboGreenwoodwidth
      !%Type float
      !%Default 0.0
      !%Section Linear Response::Kubo Greenwood
      !%Description
      !% The width applied to the Kubo Greenwood conductivity.
      !% In units of energy. Cannot be negative.
      !%End

      call parse_variable('KuboGreenwoodWidth', M_ZERO, width, units_inp%energy)
      if(width < -M_EPSILON) then
        message(1) = "KuboGreenwoodWidth cannot be negative."
        call messages_fatal(1)
      end if

      !%Variable KuboGreenwoodMaxEnergy
      !%Type float
      !%Default 2.0
      !%Section Linear Response::Kubo Greenwood
      !%Description
      !% Maximum energy for Kubo Greenwood
      !% In units of energy. 
      !%End

      call parse_variable('KuboGreenwoodMaxEnergy', CNST(2.0),maxfreq , units_inp%energy)



      !%Variable KuboGreenwoodMinEnergy
      !%Type float
      !%Default 0.001
      !%Section Linear Response::Kubo Greenwood
      !%Description
      !% Minimum energy for Kubo Greenwood
      !% In units of energy. 
      !%End

      call parse_variable('KuboGreenwoodMinEnergy', CNST(0.001),minfreq , units_inp%energy)

      !%Variable KuboGreenwoodEnergySpacing
      !%Type float
      !%Default 0.01
      !%Section Linear Response::Kubo Greenwood
      !%Description
      !% Spacing of frequencies for Kubo Greenwood
      !% In units of energy. 
      !%End

      call parse_variable('KuboGreenwoodEnergySpacing', CNST(0.01), dfreq, units_inp%energy)

      
      call restart_init(gs_restart, RESTART_GS, RESTART_TYPE_LOAD, sys%mc, ierr, mesh = sys%gr%mesh, exact = .true.)
    if(ierr == 0) then
      call states_look_and_load(gs_restart, sys%st, sys%gr)
      call restart_end(gs_restart)
    else
      call messages_write("Cannot find occupied and unoccupied states, a previous ground state calculation is required.")
      call messages_fatal()
    end if

    ASSERT(.not. sys%st%parallel_in_states)    

    mesh => sys%gr%mesh
    
    SAFE_ALLOCATE(psii(1:mesh%np_part, 1:sys%st%d%dim))
    SAFE_ALLOCATE(psij(1:mesh%np_part, 1:sys%st%d%dim))
    SAFE_ALLOCATE(gpsii(1:mesh%np, 1:sys%gr%sb%dim, 1:sys%st%d%dim))
    SAFE_ALLOCATE(gpsij(1:mesh%np, 1:sys%gr%sb%dim, 1:sys%st%d%dim))
    nfreq = nint(maxfreq / (dfreq)) + 1
    SAFE_ALLOCATE(tensor(1:mesh%sb%dim, 1:mesh%sb%dim,1:nfreq))
    SAFE_ALLOCATE(therm_tensor(1:mesh%sb%dim, 1:mesh%sb%dim,1:nfreq))
    SAFE_ALLOCATE(k11(1:mesh%sb%dim, 1:mesh%sb%dim, 1:nfreq))
    SAFE_ALLOCATE(k12(1:mesh%sb%dim, 1:mesh%sb%dim, 1:nfreq))
    SAFE_ALLOCATE(k21(1:mesh%sb%dim, 1:mesh%sb%dim, 1:nfreq))
    SAFE_ALLOCATE(k22(1:mesh%sb%dim, 1:mesh%sb%dim, 1:nfreq))
    tensor = CNST(0.0)
    k12 = CNST(0.0)
    k11 = CNST(0.0)
    k21 = CNST(0.0)
    k22 = CNST(0.0)
    therm_tensor = CNST(0.0)
    sum_occ = CNST(0.0)
    sum_cond1 = CNST(0.0)
    sum_cond2 = CNST(0.0)
    sum_cond = CNST(0.0)
    sum_k = CNST(0.0)
 

    call smear_find_fermi_energy(sys%st%smear, sys%st%eigenval, sys%st%occ, sys%st%qtot, &
      sys%st%d%nik, sys%st%nst, sys%st%d%kweights)

 !!   write(*,*)"smear", sys%st%smear%dsmear
    
    e_fermi = sys%st%smear%e_fermi
    
    do iqn = sys%st%d%kpt%start, sys%st%d%kpt%end
       
       do ist = 1, sys%st%nst
        ! get the state i and calculate the gradient
        call states_get_state(sys%st, mesh, ist, iqn, psii)
        do idim = 1, sys%st%d%dim
          call boundaries_set(sys%gr%der%boundaries, psii(:, idim))
        end do
        if(associated(hm%hm_base%phase)) then 
          call states_set_phase(sys%st%d, psii, hm%hm_base%phase(:, iqn), mesh%np_part, .false.)
       end if
        do idim = 1, sys%st%d%dim 
          call zderivatives_grad(sys%gr%der, psii(:, idim), gpsii(:, :, idim), set_bc = .false.)
        end do

         occi = sys%st%occ(ist,iqn)
         sum_occ = sum_occ + occi*sys%st%d%kweights(iqn)
        
        do jst = 1, sys%st%nst
          ! get the state j and calculate the gradient
          call states_get_state(sys%st, mesh, jst, iqn, psij)
          do idim = 1, sys%st%d%dim
            call boundaries_set(sys%gr%der%boundaries, psij(:, idim))
          end do
          if(associated(hm%hm_base%phase)) then 
            call states_set_phase(sys%st%d, psij, hm%hm_base%phase(:, iqn), mesh%np_part, .false.)
          end if
          do idim = 1, sys%st%d%dim
            call zderivatives_grad(sys%gr%der, psij(:, idim), gpsij(:, :, idim), set_bc = .false.)
          end do
          eigi = sys%st%eigenval(ist,iqn)
          eigj = sys%st%eigenval(jst,iqn)
          occi = sys%st%occ(ist,iqn)
          occj = sys%st%occ(jst,iqn)
          if (abs(eigi-eigj) < CNST(1.0e-10)) then
             df = smear_step_function_der(sys%st%smear,eigi)
             !df = step_der(sys%st%smear,eigi)
          else
             df = (occj - occi)/(eigj - eigi)
          endif

          do idir = 1, mesh%sb%dim
             do jdir = 1, mesh%sb%dim
                prod =  zmf_dotp(mesh, sys%st%d%dim, psii, gpsij(:, idir, :)) * zmf_dotp(mesh, sys%st%d%dim, psij, gpsii(:, jdir, :))
                do ifreq = 1, nfreq
                   if (abs(eigi-eigj)<CNST(1.0e-10)) then 
                      tensor(idir, jdir,ifreq) = tensor(idir, jdir,ifreq) + sys%st%d%kweights(iqn)*(CNST(-1.0)/mesh%sb%rcell_volume)* &
                           df*real(prod,REAL_PRECISION) * (CNST(0.5)*width + M_ZI*((minfreq + (ifreq-1)*dfreq)))/ (((minfreq + (ifreq-1)*dfreq))**2 + width**2/CNST(4.0)) 
                      k12(idir,jdir,ifreq) = k12(idir,jdir,ifreq) + (eigi - e_fermi)*sys%st%d%kweights(iqn)*(CNST(-1.0)/mesh%sb%rcell_volume)* &
                           df*real(prod,REAL_PRECISION) * (CNST(0.5)*width - M_ZI*(eigi-eigj - (minfreq + (ifreq-1)*dfreq)))/ ((eigi-eigj - (minfreq + (ifreq-1)*dfreq))**2 + width**2/CNST(4.0))
                      k21(idir,jdir,ifreq) = k21(idir,jdir,ifreq) + (eigj - e_fermi)*sys%st%d%kweights(iqn)*(CNST(-1.0)/mesh%sb%rcell_volume)* &
                           df*real(prod,REAL_PRECISION) * (CNST(0.5)*width - M_ZI*(eigi-eigj - (minfreq + (ifreq-1)*dfreq)))/ ((eigi-eigj - (minfreq + (ifreq-1)*dfreq))**2 + width**2/CNST(4.0))
                      k22(idir,jdir,ifreq) = k22(idir,jdir,ifreq) - (eigi - e_fermi)*(eigj - e_fermi)*sys%st%d%kweights(iqn)*(CNST(-1.0)/mesh%sb%rcell_volume)* &
                          df*real(prod,REAL_PRECISION) * (CNST(0.5)*width - M_ZI*(eigi-eigj - (minfreq + (ifreq-1)*dfreq)))/ ((eigi-eigj - (minfreq + (ifreq-1)*dfreq))**2 + width**2/CNST(4.0))
                   else 
                      tensor(idir, jdir,ifreq) = tensor(idir, jdir,ifreq) - sys%st%d%kweights(iqn)*(CNST(-1.0)/mesh%sb%rcell_volume)* &
                          df*real(prod,REAL_PRECISION) * (CNST(0.5)*width - M_ZI*(eigi-eigj - (minfreq + (ifreq-1)*dfreq)))/ ((eigi-eigj - (minfreq + (ifreq-1)*dfreq))**2 + width**2/CNST(4.0))
                      k12(idir,jdir,ifreq) = k12(idir,jdir,ifreq) + (eigi - e_fermi)*sys%st%d%kweights(iqn)*(CNST(-1.0)/mesh%sb%rcell_volume)* &
                           df*real(prod,REAL_PRECISION) * (CNST(0.5)*width - M_ZI*(eigi-eigj - (minfreq + (ifreq-1)*dfreq)))/ ((eigi-eigj - (minfreq + (ifreq-1)*dfreq))**2 + width**2/CNST(4.0))
                      k21(idir,jdir,ifreq) = k21(idir,jdir,ifreq) + (eigj - e_fermi)*sys%st%d%kweights(iqn)*(CNST(-1.0)/mesh%sb%rcell_volume)* &
                           df*real(prod,REAL_PRECISION) * (CNST(0.5)*width - M_ZI*(eigi-eigj - (minfreq + (ifreq-1)*dfreq)))/ ((eigi-eigj - (minfreq + (ifreq-1)*dfreq))**2 + width**2/CNST(4.0))
                      k22(idir,jdir,ifreq) = k22(idir,jdir,ifreq) - (eigi - e_fermi)*(eigj - e_fermi)*sys%st%d%kweights(iqn)*(CNST(-1.0)/mesh%sb%rcell_volume)* &
                          df*real(prod,REAL_PRECISION) * (CNST(0.5)*width - M_ZI*(eigi-eigj - (minfreq + (ifreq-1)*dfreq)))/ ((eigi-eigj - (minfreq + (ifreq-1)*dfreq))**2 + width**2/CNST(4.0))
                  endif
                end do
             end do !loop over jdir
          end do !loop over idir
       end do !loop over states j
    end do !loop over states i

    sum_k = sum_k + sys%st%d%kweights(iqn)
   end do !kpt loop

    !begin sum rule

   do ifreq=1,nfreq
      sum_cond1 = sum_cond1 + (tensor(1,1,ifreq) + tensor(2,2,ifreq) + tensor(3,3,ifreq))*dfreq/CNST(3.0)
      !! thermal conductivity 
      do idir = 1, mesh%sb%dim
         do jdir = 1, mesh%sb%dim
            if (abs(real(tensor(idir,jdir,ifreq),REAL_PRECISION)) .lt. CNST(1.0e-10)) then
               therm_tensor(idir,jdir,ifreq) = therm_tensor(idir,jdir,ifreq) + k22(idir,jdir,ifreq) 
            else
               therm_tensor(idir,jdir,ifreq) = therm_tensor(idir,jdir,ifreq) + (k22(idir,jdir,ifreq) - (k12(idir,jdir,ifreq)*k21(idir,jdir,ifreq)/tensor(idir,jdir,ifreq))) 
            endif
         enddo
      enddo
      
   enddo
   sum_cond2 = sum_cond2 + CNST(2.0) * mesh%sb%rcell_volume * sum_cond1 / (sum_occ * M_PI)

   
   write(*,*)"sum rule for occupations::: ", sum_occ
   write(*,*)"sum rule for conductivity:: ", sum_cond2
   write(*,*)"sum of k weights ", sum_k 

   !output
   write(dirname, '(a, a)') 'kubo_greenwood' 
   call io_mkdir(trim(dirname))

   iunit = io_open(trim(dirname)//'/electrical_conductivity', action='write')

    write(unit = iunit, iostat = ierr, fmt = '(a)') &
      '###########################################################################################################################'
    write(unit = iunit, iostat = ierr, fmt = '(8a)')  '# HEADER'
    write(unit = iunit, iostat = ierr, fmt = '(a,a,a)') &
      !'#  Energy [', trim(units_abbrev(units_out%energy)), '] Conductivity [a.u.] freq ReXX ImXX ReYY ImYY ReZZ ImZZ'
      'freq ReXX ImXX ReYY ImYY ReZZ ImZZ'
    write(unit = iunit, iostat = ierr, fmt = '(a)') &
      '###########################################################################################################################'
    do ifreq = 1, nfreq
       ww = minfreq + (ifreq-1)*dfreq
       write(unit = iunit, iostat = ierr, fmt = '(7e20.10)') ww, &
            real(tensor(1,1,ifreq)), real(tensor(2,2,ifreq)), real(tensor(3,3,ifreq))
    end do

   
   call io_close(iunit)


   !output
   write(dirname, '(a, a)') 'kubo_greenwood' 
   call io_mkdir(trim(dirname))

   iunit = io_open(trim(dirname)//'/L_coefficients', action='write')

    write(unit = iunit, iostat = ierr, fmt = '(a)') &
      '###########################################################################################################################'
    write(unit = iunit, iostat = ierr, fmt = '(8a)')  '# HEADER'
    write(unit = iunit, iostat = ierr, fmt = '(a,a,a)') &
      !'#  Energy [', trim(units_abbrev(units_out%energy)), '] Conductivity [a.u.] freq ReXX ImXX ReYY ImYY ReZZ ImZZ'
      'freq L12 L21 L22'
    write(unit = iunit, iostat = ierr, fmt = '(a)') &
      '###########################################################################################################################'
    do ifreq = 1, nfreq
       ww = minfreq + (ifreq-1)*dfreq
       write(unit = iunit, iostat = ierr, fmt = '(7e20.10)') ww, &
            k12(1,1,ifreq), k21(1,1,ifreq), k22(1,1,ifreq) 
    end do
   call io_close(iunit)



   
   !output
   write(dirname, '(a, a)') 'kubo_greenwood' 
   call io_mkdir(trim(dirname))

   iunit = io_open(trim(dirname)//'/thermal_conductivity', action='write')

    write(unit = iunit, iostat = ierr, fmt = '(a)') &
      '###########################################################################################################################'
    write(unit = iunit, iostat = ierr, fmt = '(8a)')  '# HEADER'
    write(unit = iunit, iostat = ierr, fmt = '(a,a,a)') &
      !'#  Energy [', trim(units_abbrev(units_out%energy)), '] Conductivity [a.u.] freq ReXX ImXX ReYY ImYY ReZZ ImZZ'
      'freq ReXX ImXX ReYY ImYY ReZZ ImZZ'
    write(unit = iunit, iostat = ierr, fmt = '(a)') &
      '###########################################################################################################################'
    do ifreq = 1, nfreq
       ww = (ifreq-1)*dfreq
       write(unit = iunit, iostat = ierr, fmt = '(7e20.10)') ww, &
            therm_tensor(1,1,ifreq), therm_tensor(2,2,ifreq), therm_tensor(3,3,ifreq) 
    end do

   
   call io_close(iunit)
!end output
    SAFE_DEALLOCATE_A(tensor)
    SAFE_DEALLOCATE_A(psii)
    SAFE_DEALLOCATE_A(psij)
    SAFE_DEALLOCATE_A(gpsii)
    SAFE_DEALLOCATE_A(gpsij)    
    
    POP_SUB(kubo_greewood_run)
  end subroutine kubo_greenwood_run



  FLOAT function step_der(this, energy) result(stepf)
    type(smear_t), intent(in) :: this
    FLOAT,         intent(in) :: energy

    FLOAT, parameter :: maxarg = CNST(200.0)
    FLOAT :: dsmear, xx

    PUSH_SUB(smear_step_function_der)

    dsmear = max(CNST(1e-14), this%dsmear)
    
    xx = (this%e_fermi - energy) / dsmear

    !write(*,*)"dsmear",dsmear
    !write(*,*)"xx", xx
    stepf = M_ZERO
    select case(this%method)
  
    case(SMEAR_FERMI_DIRAC)
       stepf = M_ONE / (exp(-xx) + M_ONE)
       stepf = CNST(-1.0)*stepf*(stepf-CNST(1.0))/dsmear*CNST(2.0)
       !stepf = M_ONE / (CNST(2.0) + exp(-xx) + exp(xx)) / dsmear
       
     !  if (xx > maxarg) then
     !    stepf = M_ONE
     !  else if(xx > -maxarg) then
     !     stepf = M_ONE / (M_ONE + exp(-xx))
     !     stepf = stepf*stepf*exp(-xx) / dsmear
     ! end if

    !!!write(*,*)"stepf",stepf
    case default
       call messages_not_implemented('Occupation derivatives only implemented for Fermi Dirac')
  
    end select

    POP_SUB(step_der)
  end function step_der

end module kubo_greenwood_oct_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
