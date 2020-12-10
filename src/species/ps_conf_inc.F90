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

!> This routines provides, given Z and the number of valence electron the occupations of the
!> orbitals.
!> The occupations are stored in an array of size lmax, which assumes that there is only
!> value of the principal quantum number occupied for a given angular momentum
subroutine ps_guess_atomic_occupations(zval, valcharge, lmax, ispin, occ)
 FLOAT,   intent(in)     :: zval
 FLOAT,   intent(in)     :: valcharge
 integer, intent(in)     :: lmax
 integer, intent(in)     :: ispin
 FLOAT,   intent(inout)  :: occ(:,:) !< (llmax+1, 2) 

 integer :: ll
 FLOAT :: val, x

 PUSH_SUB(ps_guess_atomic_occupations)

 val = valcharge

 ASSERT(valcharge <= zval)
 if(valcharge > 2) then
   ASSERT(lmax > 0)
 end if
 if(valcharge > 8) then
   ASSERT(lmax > 1)
 end if
 if(valcharge > 18) then
   ASSERT(lmax > 2)
 end if

 select case(int(zval))
 case(1)
   call fill_s_orbs(val, 1) ! H 1s^1
 case(2)
   call fill_s_orbs(val, 2) ! He 2s^2
 case(3) 
   call fill_s_orbs(val, 1) ! Li 1s^2 2s^1
 case(4) 
   call fill_s_orbs(val, 2) ! Be 1s^2 ; 2s^2
 case(5, 13)
   if(val > 2) call fill_s_orbs(val, 2)
   call fill_p_orbs(val, 1) ! B  1s^2 ; 2s^2 2p^1
                           ! Al      1s^2 2s^2 2p^6 ; 3s^2 3p^1 
 case(6, 14)
   if(val > 2) call fill_s_orbs(val, 2)
   call fill_p_orbs(val, 2)! C       1s^2 ; 2s^2 2p^2
                          ! Si      1s^2 2s^2 2p^6 ; 3s^2 3p^2
 case(7, 15)
   if(val > 2) call fill_s_orbs(val, 2)
   call fill_p_orbs(val, 3)! N       1s^2 ; 2s^2 2p^3
                          ! P       1s^2 2s^2 2p^6 ; 3s^2 3p^3
 case(8, 16)
   if(val > 2) call fill_s_orbs(val, 2)
   call fill_p_orbs(val, 4) ! O  1s^2 ; 2s^2 2p^4
                           ! S       1s^2 2s^2 2p^6 ; 3s^2 3p^4
 case(9, 17)
   if(val > 2) call fill_s_orbs(val, 2)
   call fill_p_orbs(val, 5)! F       1s^2 ; 2s^2 2p^5
                          ! Cl      1s^2 2s^2 2p^6 ; 3s^2 3p^5
 case(10, 18)
   if(val > 2) call fill_s_orbs(val, 2)
   call fill_p_orbs(val, 6)! Ne       1s^2 ; 2s^2 2p^6
                          ! Ar      1s^2 2s^2 2p^6 ; 3s^2 3p^6
 case(11, 19)
   if(val > 1) call fill_p_orbs(val, 6)
   call fill_s_orbs(val, 1) ! Na      1s^2 2s^2 2p^6 ; 3s^1
                           ! K       1s^2 2s^2 2p^6 3s^2 3p^6 ; 4s^1
 case(12, 20)
   if(val > 2) call fill_p_orbs(val, 6)
   call fill_s_orbs(val, 2) ! Mg      1s^2 2s^2 2p^6 ; 3s^2
                           ! Ca      1s^2 2s^2 2p^6 3s^2 3p^6 ; 4s^2
 case(21)
   if(val > 6) call fill_p_orbs(val, 6)
   if(val > 2) call fill_s_orbs(val, 2)
   call fill_d_orbs(val, 1) ! Sc      1s^2 2s^2 2p^6 3s^2 3p^6 ; 4s^2 3d^1
 case(22)
   if(val > 6) call fill_p_orbs(val, 6)
   if(val > 2) call fill_s_orbs(val, 2)
   call fill_d_orbs(val, 2) ! Ti      1s^2 2s^2 2p^6 3s^2 3p^6 ; 4s^2 3d^2
 case(23)
   if(val > 6) call fill_p_orbs(val, 6)
   if(val > 2) call fill_s_orbs(val, 2)
   call fill_d_orbs(val, 3) ! V       1s^2 2s^2 2p^6 3s^2 3p^6 ; 4s^2 3d^3
 case(24)
   if(val > 6) call fill_p_orbs(val, 6)
   if(val > 2) call fill_s_orbs(val, 2)
   call fill_d_orbs(val, 4) ! Cr      1s^2 2s^2 2p^6 3s^2 3p^6 ; 4s^2 3d^4
 case(25)
   if(val > 6) call fill_p_orbs(val, 6)
   if(val > 2) call fill_s_orbs(val, 2)
   call fill_d_orbs(val, 5) ! Mn      1s^2 2s^2 2p^6 3s^2 3p^6 ; 4s^2 3d^5
 case(26)
   if(val > 6) call fill_p_orbs(val, 6)
   if(val > 2) call fill_s_orbs(val, 2)
   call fill_d_orbs(val, 6) ! Fe      1s^2 2s^2 2p^6 3s^2 3p^6 ; 4s^2 3d^6
 case(27)
   if(val > 6) call fill_p_orbs(val, 6)
   if(val > 2) call fill_s_orbs(val, 2)
   call fill_d_orbs(val, 7) ! Co      1s^2 2s^2 2p^6 3s^2 3p^6 ; 4s^2 3d^7
 case(28)
   if(val > 6) call fill_p_orbs(val, 6)
   if(val > 2) call fill_s_orbs(val, 2)
   call fill_d_orbs(val, 8) ! Ni      1s^2 2s^2 2p^6 3s^2 3p^6 ; 4s^2 3d^8
 case(29)
   if(val > 16) then ! Cu      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 ; 4s^1
     call fill_p_orbs(val, 6)
     call fill_d_orbs(val, 10)
     call fill_s_orbs(val, 1)
   else if(val > 10) then
     call fill_d_orbs(val, 10)
     call fill_s_orbs(val, 1)
   else
     call fill_s_orbs(val, 1)
   end if
 case(30)
   if(val > 16) then ! Zn      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 ; 4s^2
     call fill_p_orbs(val, 6)
     call fill_d_orbs(val, 10)
     call fill_s_orbs(val, 2)
   else if(val > 10) then
     call fill_d_orbs(val, 10)
     call fill_s_orbs(val, 2)
   else
     call fill_s_orbs(val, 2)
   end if
 case(31)
   if(val > 10) call fill_d_orbs(val, 10)
   if(val > 2) call fill_s_orbs(val, 2)
   call fill_p_orbs(val, 1) ! Ga      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 ; 4s^2 4p^1
 case(32)
   if(val > 10) call fill_d_orbs(val, 10)
   if(val > 2) call fill_s_orbs(val, 2)
   call fill_p_orbs(val, 2) ! Ge      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 ; 4s^2 4p^2
 case(33)
   if(val > 10) call fill_d_orbs(val, 10)
   if(val > 2) call fill_s_orbs(val, 2)
   call fill_p_orbs(val, 3) ! As      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 ; 4s^2 4p^3
 case(34)
   if(val > 10) call fill_d_orbs(val, 10)
   if(val > 2) call fill_s_orbs(val, 2)
   call fill_p_orbs(val, 4) ! Se      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 ; 4s^2 4p^4
 case(35)
   if(val > 10) call fill_d_orbs(val, 10)
   if(val > 2) call fill_s_orbs(val, 2)
   call fill_p_orbs(val, 5) ! Br      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 ; 4s^2 4p^5
 case(36)
   if(val > 10) call fill_d_orbs(val, 10)
   if(val > 2) call fill_s_orbs(val, 2)
   call fill_p_orbs(val, 6) ! Kr      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 ; 4s^2 4p^6
 case(37) 
   if(val > 6) call fill_p_orbs(val, 6)
   call fill_s_orbs(val, 1) ! Rb      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 ; 5s^1
 case(38)
   if(val > 6) call fill_p_orbs(val, 6)
   call fill_s_orbs(val, 1) ! Sr      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 ; 5s^2
 case(39)
   if(val > 6) call fill_p_orbs(val, 6)
   if(val > 2) call fill_d_orbs(val, 1)
   call fill_s_orbs(val, 2) ! Y       1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 ; 4d^1 5s^2
 case(40)
   if(val > 6) call fill_p_orbs(val, 6)
   if(val > 2) call fill_d_orbs(val, 2)
   call fill_s_orbs(val, 2) ! Zr      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 ; 4d^2 5s^2
 case(41)
   if(val > 6) call fill_p_orbs(val, 6)
   if(val > 1) call fill_d_orbs(val, 4)
   call fill_s_orbs(val, 1) ! Nb      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 ; 4d^4 5s^1
 case(42)
   if(val > 6) call fill_p_orbs(val, 6)
   if(val > 1) call fill_d_orbs(val, 5)
   call fill_s_orbs(val, 1) ! Mo      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 ; 4d^5 5s^1
 case(43)
   if(val > 7) call fill_p_orbs(val, 6)
   if(val > 2) call fill_d_orbs(val, 5)
   call fill_s_orbs(val, 2) ! Tc      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 ; 4d^5 5s^2
 case(44)
   if(val > 8) call fill_p_orbs(val, 6)
   if(val > 1) call fill_d_orbs(val, 7)
   call fill_s_orbs(val, 1) ! Ru      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 ; 4d^7 5s^1
 case(45) 
   if(val > 9) call fill_p_orbs(val, 6)
   if(val > 1) call fill_d_orbs(val, 8)
   call fill_s_orbs(val, 1) ! Rh      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 ; 4d^8 5s^1
 case(46)
   if(val > 10) call fill_p_orbs(val, 6)
   call fill_d_orbs(val, 10)
 case(47)
   if(val > 16) call fill_p_orbs(val, 6)
   if(val > 10) call fill_d_orbs(val, 10)
   call fill_s_orbs(val, 1) ! Ag      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 ; 5s^1
 case(48)
   if(val > 16) call fill_p_orbs(val, 6)
   if(val > 10) call fill_d_orbs(val, 10)
   call fill_s_orbs(val, 1) ! Cd      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 ; 5s^2
 case(49)
   if(val > 10) call fill_d_orbs(val, 10)
   if(val > 1) call fill_s_orbs(val, 2)
   call fill_p_orbs(val, 1) ! In      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 ; 5s^2 5p^1
 case(50)
   if(val > 10) call fill_d_orbs(val, 10)
   if(val > 2) call fill_s_orbs(val, 2)
   call fill_p_orbs(val, 2) ! Sn      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 ; 5s^2 5p^2
 case(51)
   if(val > 10) call fill_d_orbs(val, 10)
   if(val > 3) call fill_s_orbs(val, 2)
   call fill_p_orbs(val, 3) ! Sb      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 ; 5s^2 5p^3
 case(52)
   if(val > 10) call fill_d_orbs(val, 10)
   if(val > 4) call fill_s_orbs(val, 2)
   call fill_p_orbs(val, 4) ! Te      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 ; 5s^2 5p^4
 case(53)
   if(val > 10) call fill_d_orbs(val, 10)
   if(val > 5) call fill_s_orbs(val, 2)
   call fill_p_orbs(val, 5) ! I       1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 ; 5s^2 5p^5
 case(54)
   if(val > 10) call fill_d_orbs(val, 10)
   if(val > 6) call fill_s_orbs(val, 2)
   call fill_p_orbs(val, 6) ! Xe      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 ; 5s^2 5p^6
 case(55)
   if(val > 6) call fill_p_orbs(val, 6)
   call fill_s_orbs(val, 1) ! Cs      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 5s^2 5p^6 ; 6s^1
 case(56)
   if(val > 6) call fill_p_orbs(val, 6)
   call fill_s_orbs(val, 2) ! Ba      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 5s^2 5p^6 ; 6s^2
 case(57)
   if(val > 6) call fill_p_orbs(val, 6)
   if(val > 2) call fill_d_orbs(val, 1)
   call fill_s_orbs(val, 2) ! La      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 5s^2 5p^6 ; 5d^1 6s^2
 case(58)
   if(val > 6) call fill_p_orbs(val, 6)
   if(val > 3) call fill_f_orbs(val, 1)
   if(val > 2) call fill_d_orbs(val, 1)
   call fill_s_orbs(val, 2) ! Ce      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 5s^2 5p^6 ; 4f^1 5d^1 6s^2
 case(59)
   if(val > 6) call fill_p_orbs(val, 6)
   if(val > 2) call fill_f_orbs(val, 3)
   call fill_s_orbs(val, 2) ! Pr      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 5s^2 5p^6 ; 4f^3 6s^2
 case(60)
   if(val > 6) call fill_p_orbs(val, 6)
   if(val > 2) call fill_f_orbs(val, 4)
   call fill_s_orbs(val, 2) ! Nd      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 5s^2 5p^6 ; 4f^4 6s^2
 case(61)
   if(val > 7) call fill_p_orbs(val, 6)
   if(val > 2) call fill_f_orbs(val, 5)
   call fill_s_orbs(val, 2) ! Pm      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 5s^2 5p^6 ; 4f^5 6s^2
 case(62)
   if(val > 8) call fill_p_orbs(val, 6)
   if(val > 2) call fill_f_orbs(val, 6)
   call fill_s_orbs(val, 2) ! Sm      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 5s^2 5p^6 ; 4f^6 6s^2
 case(63)
   if(val > 9) call fill_p_orbs(val, 6)
   if(val > 2) call fill_f_orbs(val, 7)
   call fill_s_orbs(val, 2) ! Eu      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 5s^2 5p^6 ; 4f^7 6s^2
 case(64)
   if(val > 10) call fill_p_orbs(val, 6)
   if(val > 3) call fill_f_orbs(val, 7)
   if(val > 2) call fill_d_orbs(val, 1)
   call fill_s_orbs(val, 2) ! Gd      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 5s^2 5p^6 ; 4f^7 5d^1 6s^2
 case(65)
   if(val > 11) call fill_p_orbs(val, 6)
   if(val > 2) call fill_f_orbs(val, 9)
   call fill_s_orbs(val, 2) ! Tb      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 5s^2 5p^6 ; 4f^9 6s^2
 case(66)
   if(val > 12) call fill_p_orbs(val, 6)
   if(val > 2) call fill_f_orbs(val, 10)
   call fill_s_orbs(val, 2) ! Dy      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 5s^2 5p^6 ; 4f^10 6s^2
 case(67)
   if(val > 13) call fill_p_orbs(val, 6)
   if(val > 2) call fill_f_orbs(val, 11)
   call fill_s_orbs(val, 2) ! Ho      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 5s^2 5p^6 ; 4f^11 6s^2
 case(68)
   if(val > 14) call fill_p_orbs(val, 6)
   if(val > 2) call fill_f_orbs(val, 12)
   call fill_s_orbs(val, 2) ! Er      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 5s^2 5p^6 ; 4f^12 6s^2
 case(69)
   if(val > 15) call fill_p_orbs(val, 6)
   if(val > 2) call fill_f_orbs(val, 13)
   call fill_s_orbs(val, 2) ! Tm      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 5s^2 5p^6 ; 4f^13 6s^2
 case(70)
   if(val > 16) call fill_p_orbs(val, 6)
   if(val > 2) call fill_f_orbs(val, 14)
   call fill_s_orbs(val, 2) ! Yb      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 5s^2 5p^6 ; 4f^14 6s^2
 case(71)
   if(val > 20) call fill_p_orbs(val, 6)
   if(val > 14) call fill_f_orbs(val, 14)
   if(val > 2) call fill_d_orbs(val, 1)
   call fill_s_orbs(val, 2) ! Lu      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 5s^2 5p^6 4f^14 ; 5d^1 6s^2
 case(72)
   if(val > 20) call fill_p_orbs(val, 6)
   if(val > 14) call fill_f_orbs(val, 14)
   if(val > 2) call fill_d_orbs(val, 2)
   call fill_s_orbs(val, 2) ! Hf      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 5s^2 5p^6 4f^14 ; 5d^2 6s^2
 case(73)
   if(val > 20) call fill_p_orbs(val, 6)
   if(val > 14) call fill_f_orbs(val, 14)
   if(val > 2) call fill_d_orbs(val, 3)
   call fill_s_orbs(val, 2) ! Ta      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 5s^2 5p^6 4f^14 ; 5d^3 6s^2
 case(74)
   if(val > 20) call fill_p_orbs(val, 6)
   if(val > 14) call fill_f_orbs(val, 14)
   if(val > 2) call fill_d_orbs(val, 4)
   call fill_s_orbs(val, 2) ! W       1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 5s^2 5p^6 4f^14 ; 5d^4 6s^2
 case(75)
   if(val > 21) call fill_p_orbs(val, 6)
   if(val > 14) call fill_f_orbs(val, 14)
   if(val > 2) call fill_d_orbs(val, 5)
   call fill_s_orbs(val, 2) ! Re      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 5s^2 5p^6 4f^14 ; 5d^5 6s^2
 case(76)
   if(val > 22) call fill_p_orbs(val, 6)
   if(val > 14) call fill_f_orbs(val, 14)
   if(val > 2) call fill_d_orbs(val, 6)
   call fill_s_orbs(val, 2) ! Os      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 5s^2 5p^6 4f^14 ; 5d^6 6s^2
 case(77)
   if(val > 23) call fill_p_orbs(val, 6)
   if(val > 14) call fill_f_orbs(val, 14)
   if(val > 2) call fill_d_orbs(val, 7)
   call fill_s_orbs(val, 2) ! Ir      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 5s^2 5p^6 4f^14 ; 5d^7 6s^2
 case(78)
   if(val > 24) call fill_p_orbs(val, 6)
   if(val > 14) call fill_f_orbs(val, 14)
   if(val > 1) call fill_d_orbs(val, 9)
   call fill_s_orbs(val, 1) ! Pt      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 5s^2 5p^6 4f^14 ; 5d^9 6s^1
 case(79)
   if(val > 25) call fill_p_orbs(val, 6)
   if(val > 14) call fill_f_orbs(val, 14)
   if(val > 1) call fill_d_orbs(val, 10)
   call fill_s_orbs(val, 1) ! Au      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 5s^2 5p^6 4f^14 5d^10 ; 6s^1
 case(80)
   if(val > 25) call fill_p_orbs(val, 6)
   if(val > 14) call fill_f_orbs(val, 14)
   if(val > 2) call fill_d_orbs(val, 10)
   call fill_s_orbs(val, 2) ! Hg      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 5s^2 5p^6 4f^14 5d^10 ; 6s^2
 case(81)
   if(val > 24) call fill_f_orbs(val, 14)
   if(val > 10) call fill_d_orbs(val, 10)
   if(val > 1) call fill_s_orbs(val, 2)
   call fill_p_orbs(val, 1) ! Tl      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 5s^2 5p^6 4f^14 5d^10 ; 6s^2 6p^1
 case(82)
   if(val > 24) call fill_f_orbs(val, 14)
   if(val > 10) call fill_d_orbs(val, 10)
   if(val > 2) call fill_s_orbs(val, 2)
   call fill_p_orbs(val, 2) ! Pb      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 5s^2 5p^6 4f^14 5d^10 ; 6s^2 6p^2
 case(83)
   if(val > 24) call fill_f_orbs(val, 14)
   if(val > 10) call fill_d_orbs(val, 10)
   if(val > 3) call fill_s_orbs(val, 2)
   call fill_p_orbs(val, 3) ! Bi      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 5s^2 5p^6 4f^14 5d^10 ; 6s^2 6p^3
 case(84)
   if(val > 24) call fill_f_orbs(val, 14)
   if(val > 10) call fill_d_orbs(val, 10)
   if(val > 4) call fill_s_orbs(val, 2)
   call fill_p_orbs(val, 4) ! Po      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 5s^2 5p^6 4f^14 5d^10 ; 6s^2 6p^4
 case(85)
   if(val > 24) call fill_f_orbs(val, 14)
   if(val > 10) call fill_d_orbs(val, 10)
   if(val > 5) call fill_s_orbs(val, 2)
   call fill_p_orbs(val, 5) ! At      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 5s^2 5p^6 4f^14 5d^10 ; 6s^2 6p^5
 case(86)
   if(val > 24) call fill_f_orbs(val, 14)
   if(val > 10) call fill_d_orbs(val, 10)
   if(val > 6) call fill_s_orbs(val, 2)
   call fill_p_orbs(val, 6) ! Rn      1s^2 2s^2 2p^6 3s^2 3p^6 3d^10 4s^2 4p^6 4d^10 5s^2 5p^6 4f^14 5d^10 ; 6s^2 6p^6

  end select

  !If we attributed all the electrons, everything went fine
  if(val < M_EPSILON ) then
    !In case of spin-polarized calculations, we properly distribute the electrons
    if(ispin == 2) then
      do ll = 0, lmax
        x = occ(ll+1, 1)
        occ(ll+1, 1) = min(x, TOFLOAT(2*ll+1))
        occ(ll+1, 2) = x - occ(ll+1,1) 
      end do
    end if
  else
    occ = M_ZERO
  end if

  POP_SUB(ps_guess_atomic_occupations)

  contains
    subroutine fill_s_orbs(val,max_occ)
      FLOAT, intent(inout) :: val
      integer, intent(in)    :: max_occ

      occ(1,1) = min(val, real(max_occ))
      val = val - occ(1,1)
    end subroutine fill_s_orbs
 
    subroutine fill_p_orbs(val,max_occ)
      FLOAT, intent(inout) :: val
      integer, intent(in)    :: max_occ

      occ(2,1) = min(val, real(max_occ))
      val = val - occ(2,1)
    end subroutine fill_p_orbs

    subroutine fill_d_orbs(val,max_occ)
      FLOAT, intent(inout) :: val
      integer, intent(in)    :: max_occ

      occ(3,1) = min(val, real(max_occ))
      val = val - occ(3,1)
    end subroutine fill_d_orbs
    
    subroutine fill_f_orbs(val,max_occ)
      FLOAT, intent(inout) :: val
      integer, intent(in)    :: max_occ

      occ(4,1) = min(val, real(max_occ))
      val = val - occ(4,1)
    end subroutine fill_f_orbs

end subroutine ps_guess_atomic_occupations
