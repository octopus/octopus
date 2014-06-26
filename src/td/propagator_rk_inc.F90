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
!! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
!! 02110-1301, USA.
!!
!! $Id$


    subroutine td_runge_kutta2()
#ifndef HAVE_SPARSKIT
      ASSERT(.false.)
#else
      integer :: np_part, np, kp1, kp2, st1, st2, nspin, ik, ist, idim, j, ip
      integer :: i
      FLOAT :: dres
      CMPLX, allocatable :: zphi(:, :, :, :)
      CMPLX, allocatable :: zpsi(:), rhs(:)
      CMPLX, allocatable :: k2(:, :, :, :), oldk2(:, :, :, :), rhs1(:, :, :, :)
      CMPLX, allocatable :: inhpsi(:)
      type(ion_state_t) :: ions_state

      PUSH_SUB(propagator_dt.td_runge_kutta2)

      st1 = st%st_start
      st2 = st%st_end
      kp1 = st%d%kpt%start
      kp2 = st%d%kpt%end
      np_part = gr%mesh%np_part
      np = gr%mesh%np
      nspin = hm%d%nspin
      move_ions_op = ion_dynamics_ions_move(ions)

      sp_np = np
      sp_dim = st%d%dim
      sp_st1 = st1
      sp_st2 = st2
      sp_kp1 = kp1
      sp_kp2 = kp2
      sp_parallel = st%parallel_in_states .or. st%d%kpt%parallel
      if(sp_parallel) sp_comm = st%st_kpt_mpi_grp%comm

      ! define pointer and variables for usage in td_rk2op, td_rk2opt routines
      grid_p    => gr
      hm_p      => hm
      tr_p      => tr
      st_p      => st
      dt_op = dt
      t_op  = time - dt/M_TWO
      dim_op = st%d%dim
      xc_p      => ks%xc

      SAFE_ALLOCATE(k2(1:np_part, 1:st%d%dim, st1:st2, kp1:kp2))
      SAFE_ALLOCATE(oldk2(1:np_part, 1:st%d%dim, st1:st2, kp1:kp2))
      SAFE_ALLOCATE(zphi(1:gr%mesh%np_part, st%d%dim, st1:st2, kp1:kp2))
      SAFE_ALLOCATE(rhs1(1:np_part, 1:st%d%dim, st1:st2, kp1:kp2))
      SAFE_ALLOCATE(rhs(1:tr%tdsk_size))
      SAFE_ALLOCATE(zpsi(1:tr%tdsk_size))
      SAFE_ALLOCATE(vhxc1_op(np, nspin))

      ! First, we get the state that we want to propagate. For the moment being, only one state.
      do ik = kp1, kp2
        do ist = st1, st2
          call states_get_state(st, gr%mesh, ist, ik, zphi(:, :, ist, ik))
        end do
      end do

      if(hamiltonian_oct_exchange(hm)) call hamiltonian_prepare_oct_exchange(hm, gr%mesh, zphi, ks%xc)
      call hamiltonian_update(hm, gr%mesh, time = time-dt)
      rhs1 = M_z0
      do ik = kp1, kp2
        do ist = st1, st2
          call zhamiltonian_apply(hm_p, grid_p%der, zphi(:, :, ist, ik), rhs1(:, :, ist, ik), ist, ik, time -dt)
          if(hamiltonian_inh_term(hm)) then
            SAFE_ALLOCATE(inhpsi(1:gr%mesh%np))
            do idim = 1, st%d%dim
              call states_get_state(hm%inh_st, gr%mesh, idim, ist, ik, inhpsi)
              forall(ip = 1:gr%mesh%np) rhs1(ip, idim, ist, ik) = rhs1(ip, idim, ist, ik) + M_zI * inhpsi(ip)
            end do
            SAFE_DEALLOCATE_A(inhpsi)
          end if
        end do
      end do
      do ik = kp1, kp2
        do ist = st1, st2
          if(hamiltonian_oct_exchange(hm)) then
            call zoct_exchange_operator(hm, gr%der, rhs1(:, :, ist, ik), ist, ik)
          end if
        end do
      end do

      rhs1 = zphi - M_zI * M_HALF * dt * rhs1
      k2 = zphi

      i = 1
      do
        oldk2 = k2

        ! Set the Hamiltonian at the final time of the propagation
        if(.not.hamiltonian_oct_exchange(hm_p)) then
          do ik = kp1, kp2
            do ist = st1, st2
              call states_set_state(st, gr%mesh, ist, ik, k2(:, :, ist, ik))
            end do
          end do
          call density_calc(st, gr, st%rho)
          call v_ks_calc(ks, hm, st, geo)
        end if
        if(ion_dynamics_ions_move(ions)) then
          call ion_dynamics_save_state(ions, geo, ions_state)
          call ion_dynamics_propagate(ions, gr%sb, geo, time, dt)
          call hamiltonian_epot_generate(hm, gr, geo, st, time = time)
          vpsl1_op = hm%ep%vpsl
        end if
        call hamiltonian_update(hm, gr%mesh, time = time)
        if (i==1) then
          vhxc1_op(:, :) = tr%v_old(:, :, 0)
          i = i + 1
        else
          vhxc1_op = hm%vhxc
        end if
        t_op  = time

        if(ion_dynamics_ions_move(ions)) call ion_dynamics_restore_state(ions, geo, ions_state)

        j = 1
        do ik = kp1, kp2
          do ist = st1, st2
            do idim = 1, st%d%dim
              rhs(j:j+np-1) = rhs1(1:np, idim, ist, ik)
              j = j + np
            end do
          end do
        end do

        ! Now we populate an initial guess for the equation.
        j = 1
        do ik = kp1, kp2
          do ist = st1, st2
            do idim = 1, st%d%dim
              zpsi(j:j+np-1) = k2(1:np, idim, ist, ik)
              j = j + np
            end do
          end do
        end do

        t_op  = time - dt
        call zsparskit_solver_run(tr%tdsk, td_rk2op, td_rk2opt, zpsi, rhs)
 
        k2 = M_z0
        j = 1
        do ik = kp1, kp2
          do ist = st1, st2
            do idim = 1, st%d%dim
              k2(1:np, idim, ist, ik) = zpsi(j:j+np-1)
              j = j + np
            end do
          end do
        end do

        dres = M_ZERO
        do ik = kp1, kp2
          do ist = st1, st2
            do idim = 1, st%d%dim
              dres = dres + zmf_nrm2(gr%mesh, k2(:, idim, ist, ik) - oldk2(:, idim, ist, ik))
            end do
          end do
        end do

        if(dres < tr%scf_threshold) exit
      end do

      zphi = k2
      do ik = kp1, kp2
        do ist = st1, st2
          call states_set_state(st, gr%mesh, ist, ik, zphi(:, :, ist, ik))
        end do
      end do


      SAFE_DEALLOCATE_A(k2)
      SAFE_DEALLOCATE_A(oldk2)
      SAFE_DEALLOCATE_A(zphi)
      SAFE_DEALLOCATE_A(rhs1)
      SAFE_DEALLOCATE_A(zpsi)
      SAFE_DEALLOCATE_A(rhs)
      SAFE_DEALLOCATE_A(vhxc1_op)
      POP_SUB(propagator_dt.td_runge_kutta2)
#endif
    end subroutine td_runge_kutta2


    subroutine td_runge_kutta4()
#ifndef HAVE_SPARSKIT
      ASSERT(.false.)
#else
      logical :: converged
      integer :: np_part, np, iter, idim, ip, ist, ik, j, kp1, kp2, st1, st2, nspin
      FLOAT :: dres
      CMPLX, allocatable :: inhpsi(:)
      CMPLX, allocatable :: zpsi(:), rhs(:)
      CMPLX, allocatable :: zphi(:, :, :, :)
      type(ion_state_t) :: ions_state

      FLOAT :: a(2, 2), c(2), b(2)
      CMPLX, allocatable :: k1(:, :, :, :), k2(:, :, :, :), oldk1(:, :, :, :), &
                            oldk2(:, :, :, :), yn1(:, :, :, :), yn2(:, :, :, :), &
                            rhs1(:, :, :, :), rhs2(:, :, :, :)

      PUSH_SUB(propagator_dt.td_runge_kutta4)

      c(1) = M_HALF - sqrt(M_THREE)/M_SIX
      c(2) = M_HALF + sqrt(M_THREE)/M_SIX
      b(1) = M_HALF
      b(2) = M_HALF
      a(1, 1) = M_FOURTH
      a(1, 2) = M_FOURTH - sqrt(M_THREE)/M_SIX
      a(2, 1) = M_FOURTH + sqrt(M_THREE)/M_SIX
      a(2, 2) = M_FOURTH

      st1 = st%st_start
      st2 = st%st_end
      kp1 = st%d%kpt%start
      kp2 = st%d%kpt%end
      np_part = gr%mesh%np_part
      np = gr%mesh%np
      nspin = hm%d%nspin
      move_ions_op = ion_dynamics_ions_move(ions)

      sp_np = np
      sp_dim = st%d%dim
      sp_st1 = st1
      sp_st2 = st2
      sp_kp1 = kp1
      sp_kp2 = kp2
      sp_parallel = st%parallel_in_states .or. st%d%kpt%parallel
      if(sp_parallel) sp_comm = st%st_kpt_mpi_grp%comm

      ! define pointer and variables for usage in td_rk4op, td_rk4opt routines
      grid_p    => gr
      hm_p      => hm
      tr_p      => tr
      st_p      => st
      dt_op = dt
      t_op  = time - dt/M_TWO
      dim_op = st%d%dim

      SAFE_ALLOCATE(vhxc1_op(np, nspin))
      SAFE_ALLOCATE(vhxc2_op(np, nspin))
      SAFE_ALLOCATE(vpsl1_op(np))
      SAFE_ALLOCATE(vpsl2_op(np))
      SAFE_ALLOCATE(k1(1:np_part, 1:st%d%dim, st1:st2, kp1:kp2))
      SAFE_ALLOCATE(k2(1:np_part, 1:st%d%dim, st1:st2, kp1:kp2))
      SAFE_ALLOCATE(oldk1(1:np_part, 1:st%d%dim, st1:st2, kp1:kp2))
      SAFE_ALLOCATE(oldk2(1:np_part, 1:st%d%dim, st1:st2, kp1:kp2))
      SAFE_ALLOCATE(yn1(1:np_part, 1:st%d%dim, st1:st2, kp1:kp2))
      SAFE_ALLOCATE(yn2(1:np_part, 1:st%d%dim, st1:st2, kp1:kp2))
      SAFE_ALLOCATE(rhs1(1:np_part, 1:st%d%dim, st1:st2, kp1:kp2))
      SAFE_ALLOCATE(rhs2(1:np_part, 1:st%d%dim, st1:st2, kp1:kp2))
      SAFE_ALLOCATE(rhs(1:tr%tdsk_size))
      SAFE_ALLOCATE(zpsi(1:tr%tdsk_size))
      SAFE_ALLOCATE(zphi(1:gr%mesh%np_part, st%d%dim, st1:st2, kp1:kp2))

      ! First, we get the state that we want to propagate. For the moment being, only one state.
      do ik = kp1, kp2
        do ist = st1, st2
          call states_get_state(st, gr%mesh, ist, ik, zphi(:, :, ist, ik))
        end do
      end do
      k1 = M_z0
      k2 = M_z0

      do

        oldk1 = k1
        oldk2 = k2

        yn1 = zphi + a(1, 1) * k1 + a(1, 2) * k2
        yn2 = zphi + a(2, 1) * k1 + a(2, 2) * k2

        ! Set the Hamiltonian at time-dt + c(1) * dt
        do ik = kp1, kp2
          do ist = st1, st2
            call states_set_state(st, gr%mesh, ist, ik, yn1(:, :, ist, ik))
          end do
        end do
        call density_calc(st, gr, st%rho)
        call v_ks_calc(ks, hm, st, geo)
        if(ion_dynamics_ions_move(ions)) then
          call ion_dynamics_save_state(ions, geo, ions_state)
          call ion_dynamics_propagate(ions, gr%sb, geo, time - dt + c(1)*dt, c(1)*dt)
          call hamiltonian_epot_generate(hm, gr, geo, st, time = time - dt + c(1)*dt)
          vpsl1_op = hm%ep%vpsl
        end if
        call hamiltonian_update(hm, gr%mesh, time = time - dt + c(1)*dt)
        vhxc1_op = hm%vhxc
        t_op  = time - dt + c(1) * dt
        rhs1 = M_z0
        do ik = kp1, kp2
          do ist = st1, st2
            call zhamiltonian_apply(hm_p, grid_p%der, zphi(:, :, ist, ik), rhs1(:, :, ist, ik), ist, ik, t_op)
            if(hamiltonian_inh_term(hm)) then
              SAFE_ALLOCATE(inhpsi(1:gr%mesh%np))
              do idim = 1, st%d%dim
                call states_get_state(hm%inh_st, gr%mesh, idim, ist, ik, inhpsi)
                forall(ip = 1:gr%mesh%np) rhs1(ip, idim, ist, ik) = rhs1(ip, idim, ist, ik) + M_zI * inhpsi(ip)
              end do
              SAFE_DEALLOCATE_A(inhpsi)
            end if
          end do
        end do
        rhs1 = - M_zI * dt * rhs1
        if(ion_dynamics_ions_move(ions)) call ion_dynamics_restore_state(ions, geo, ions_state)

        ! Set the Hamiltonian at time-dt + c(2) * dt
        do ik = kp1, kp2
          do ist = st1, st2
            call states_set_state(st, gr%mesh, ist, ik, yn2(:, :, ist, ik))
          end do
        end do
        call density_calc(st, gr, st%rho)
        call v_ks_calc(ks, hm, st, geo)
        if(ion_dynamics_ions_move(ions)) then
          call ion_dynamics_save_state(ions, geo, ions_state)
          call ion_dynamics_propagate(ions, gr%sb, geo, time - dt + c(2)*dt, c(2)*dt)
          call hamiltonian_epot_generate(hm, gr, geo, st, time = time - dt + c(2)*dt)
          vpsl2_op = hm%ep%vpsl
        end if
        call hamiltonian_update(hm, gr%mesh, time = time - dt + c(2)*dt)
        vhxc2_op = hm%vhxc
        t_op  = time - dt + c(2) * dt
        rhs2 = M_z0
        do ik = kp1, kp2
          do ist = st1, st2
            call zhamiltonian_apply(hm_p, grid_p%der, zphi(:, :, ist, ik), rhs2(:, :, ist, ik), ist, ik, t_op)
            if(hamiltonian_inh_term(hm)) then
              SAFE_ALLOCATE(inhpsi(1:gr%mesh%np))
              do idim = 1, st%d%dim
                call states_get_state(hm%inh_st, gr%mesh, idim, ist, ik, inhpsi)
               forall(ip = 1:gr%mesh%np) rhs2(ip, idim, ist, ik) = rhs2(ip, idim, ist, ik) + M_zI * inhpsi(ip)
              end do
              SAFE_DEALLOCATE_A(inhpsi)
            end if
          end do
        end do
        rhs2 = -M_zI * dt * rhs2
        if(ion_dynamics_ions_move(ions)) call ion_dynamics_restore_state(ions, geo, ions_state)

        j = 1
        do ik = kp1, kp2
          do ist = st1, st2
            do idim = 1, st%d%dim
              rhs(j:j+np-1) = rhs1(1:np, idim, ist, ik)
              j = j + np
            end do
          end do
        end do
        do ik = kp1, kp2
          do ist = st1, st2
            do idim = 1, st%d%dim
              rhs(j:j+np-1) = rhs2(1:np, idim, ist, ik)
              j = j + np
            end do
          end do
        end do

        ! Now we populate an initial guess for the equation.
        j = 1
        do ik = kp1, kp2
          do ist = st1, st2
            do idim = 1, st%d%dim
              zpsi(j:j+np-1) = k1(1:np, idim, ist, ik)
              j = j + np
            end do
          end do
        end do
        do ik = kp1, kp2
          do ist = st1, st2
            do idim = 1, st%d%dim
              zpsi(j:j+np-1) = k2(1:np, idim, ist, ik)
              j = j + np
            end do
          end do
        end do

        t_op  = time - dt
        call zsparskit_solver_run(tr%tdsk, td_rk4op, td_rk4opt, zpsi, rhs)

        k1 = M_z0
        k2 = M_z0
        j = 1
        do ik = kp1, kp2
          do ist = st1, st2
            do idim = 1, st%d%dim
              k1(1:np, idim, ist, ik) = zpsi(j:j+np-1)
              j = j + np
            end do
          end do
        end do
        do ik = kp1, kp2
          do ist = st1, st2
            do idim = 1, st%d%dim
              k2(1:np, idim, ist, ik) = zpsi(j:j+np-1)
              j = j + np
            end do
          end do
        end do

        dres = M_ZERO
        do ik = kp1, kp2
          do ist = st1, st2
            do idim = 1, st%d%dim
              dres = dres + zmf_nrm2(gr%mesh, k1(:, idim, ist, ik) - oldk1(:, idim, ist, ik))
              dres = dres + zmf_nrm2(gr%mesh, k2(:, idim, ist, ik) - oldk2(:, idim, ist, ik))
            end do
          end do
        end do
        if(sp_parallel) call comm_allreduce(sp_comm, dres)
        !write(*, *) 'Residual = ', dres

        if(dres < tr%scf_threshold) exit
      end do


      zphi = zphi + b(1) * k1 + b(2) * k2
      do ik = kp1, kp2
        do ist = st1, st2
          call states_set_state(st, gr%mesh, ist, ik, zphi(:, :, ist, ik))
        end do
      end do

      SAFE_DEALLOCATE_A(rhs1)
      SAFE_DEALLOCATE_A(rhs2)
      SAFE_DEALLOCATE_A(k1)
      SAFE_DEALLOCATE_A(k2)
      SAFE_DEALLOCATE_A(oldk1)
      SAFE_DEALLOCATE_A(oldk2)
      SAFE_DEALLOCATE_A(yn1)
      SAFE_DEALLOCATE_A(yn2)
      SAFE_DEALLOCATE_A(vhxc1_op)
      SAFE_DEALLOCATE_A(vhxc2_op)
      SAFE_DEALLOCATE_A(vpsl1_op)
      SAFE_DEALLOCATE_A(vpsl2_op)
      SAFE_DEALLOCATE_A(zpsi)
      SAFE_DEALLOCATE_A(rhs)
      POP_SUB(propagator_dt.td_runge_kutta4)
#endif
    end subroutine td_runge_kutta4

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
