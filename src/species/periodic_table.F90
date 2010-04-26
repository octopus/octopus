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
!! $Id: ps.F90 2307 2006-07-29 00:50:22Z appel $

#include "global.h"

module periodic_table_m
  implicit none

  private
  public ::  &
    pt_number_from_symbol
  
  integer, public, parameter :: pt_n_elements = 103

  character(len=3), public :: pt_symbol(pt_n_elements) = (/      &
    "H  ", "He ", "Li ", "Be ", "B  ", "C  ", "N  ", "O  ",   &
    "F  ", "Ne ", "Na ", "Mg ", "Al ", "Si ", "P  ", "S  ",   &
    "Cl ", "Ar ", "K  ", "Ca ", "Sc ", "Ti ", "V  ", "Cr ",   &
    "Mn ", "Fe ", "Co ", "Ni ", "Cu ", "Zn ", "Ga ", "Ge ",   &
    "As ", "Se ", "Br ", "Kr ", "Rb ", "Sr ", "Y  ", "Zr ",   &
    "Nb ", "Mo ", "Tc ", "Ru ", "Rh ", "Pd ", "Ag ", "Cd ",   &
    "In ", "Sn ", "Sb ", "Te ", "I  ", "Xe ", "Cs ", "Ba ",   &
    "La ", "Ce ", "Pr ", "Nd ", "Pm ", "Sm ", "Eu ", "Gd ",   &
    "Tb ", "Dy ", "Ho ", "Er ", "Tm ", "Yb ", "Lu ", "Hf ",   &
    "Ta ", "W  ", "Re ", "Os ", "Ir ", "Pt ", "Au ", "Hg ",   &
    "Tl ", "Pb ", "Bi ", "Po ", "At ", "Rn ", "Fr ", "Ra ",   &
    "Ac ", "Th ", "Pa ", "U  ", "Np ", "Pu ", "Am ", "Cm ",   &
    "Bk ", "Cf ", "Es ", "Fm ", "Md ", "No ", "Lw "           &
    /)

  character(len=12), public :: pt_names(pt_n_elements) = (/             &
    "Hydrogen    ", "Helium      ", "Lithium     ", "Beryllium   ",  &
    "Boron       ", "Carbon      ", "Nitrogen    ", "Oxygen      ",  &
    "Fluorine    ", "Neon        ", "Sodium      ", "Magnesium   ",  &
    "Aluminium   ", "Silicon     ", "Phosphorus  ", "Sulfur      ",  &
    "Chlorine    ", "Argon       ", "Potassium   ", "Calcium     ",  &
    "Scandium    ", "Titanium    ", "Vanadium    ", "Chromium    ",  &
    "Manganese   ", "Iron        ", "Cobalt      ", "Nickel      ",  &
    "Copper      ", "Zinc        ", "Gallium     ", "Germanium   ",  &
    "Arsenic     ", "Selenium    ", "Bromine     ", "Krypton     ",  &
    "Rubidium    ", "Strontium   ", "Yttrium     ", "Zirconium   ",  &
    "Niobium     ", "Molybdenum  ", "Technetium  ", "Ruthenium   ",  &
    "Rhodium     ", "Palladium   ", "Silver      ", "Cadmium     ",  &
    "Indium      ", "Tin         ", "Antimony    ", "Tellurium   ",  &
    "Iodine      ", "Xenon       ", "Cesium      ", "Barium      ",  &
    "Lanthanum   ", "Cerium      ", "Praseodymium", "Neodymium   ",  &
    "Promethium  ", "Samarium    ", "Europium    ", "Gadolinium  ",  &
    "Terbium     ", "Dysprosium  ", "Holmium     ", "Erbium      ",  &
    "Thulium     ", "Ytterbium   ", "Lutetium    ", "Hafnium     ",  &
    "Tantalum    ", "Wolfram     ", "Rhenium     ", "Osmium      ",  &
    "Iridium     ", "Platinum    ", "Gold        ", "Mercury     ",  &
    "Thallium    ", "Lead        ", "Bismuth     ", "Polonium    ",  &
    "Astatine    ", "Radon       ", "Francium    ", "Radium      ",  &
    "Actinium    ", "Thorium     ", "Protactinium", "Uranium     ",  &
    "Neptunium   ", "Plutonium   ", "Americium   ", "Curium      ",  &
    "Berkelium   ", "Californium ", "Einsteinium ", "Fermium     ",  &
    "Mendelevium ", "Nobelium    ", "Lawrencium  "                   &
    /)

  FLOAT, public :: pt_atomic_mass(pt_n_elements) = (/  &
    CNST(  1.00794),   CNST(  4.002602), CNST(  6.941),    CNST(  9.01218),  &
    CNST( 10.81),      CNST( 12.011),    CNST( 14.00674),  CNST( 15.9994),   &
    CNST( 18.9984032), CNST( 20.1797),   CNST( 22.989768), CNST( 24.3050),   &
    CNST( 26.98154),   CNST( 28.0855),   CNST( 30.973762), CNST( 32.066),    &
    CNST( 35.4527),    CNST( 39.948),    CNST( 39.0983),   CNST( 40.078),    &
    CNST( 44.955910),  CNST( 47.88),     CNST( 50.9415),   CNST( 51.9961),   &
    CNST( 54.93805),   CNST( 55.847),    CNST( 58.93320),  CNST( 58.69),     &
    CNST( 63.546),     CNST( 65.39),     CNST( 69.723),    CNST( 72.61),     &
    CNST( 74.9216),    CNST( 78.96),     CNST( 79.904),    CNST( 83.80),     &
    CNST( 85.4678),    CNST( 87.62),     CNST( 88.90585),  CNST( 91.224),    &
    CNST( 92.90638),   CNST( 95.94),     CNST( 98.9062),   CNST(101.07),     &
    CNST(102.90550),   CNST(106.42),     CNST(107.8682),   CNST(112.411),    &
    CNST(114.82),      CNST(118.710),    CNST(121.75),     CNST(127.60),     &
    CNST(126.90447),   CNST(131.29),     CNST(132.9054),   CNST(137.327),    &
    CNST(138.9055),    CNST(140.115),    CNST(140.90765),  CNST(144.24),     &
    CNST(145.),        CNST(150.36),     CNST(151.965),    CNST(157.25),     &
    CNST(158.92534),   CNST(162.50),     CNST(164.93032),  CNST(167.26),     &
    CNST(168.93421),   CNST(173.04),     CNST(174.967),    CNST(178.49),     &
    CNST(180.9479),    CNST(183.85),     CNST(186.207),    CNST(190.2),      &
    CNST(192.22),      CNST(195.08),     CNST(196.96654),  CNST(200.59),     &
    CNST(204.3833),    CNST(207.2),      CNST(208.9804),   CNST(209.),       &
    CNST(210.),        CNST(222.),       CNST(223.),       CNST(226.0254),   &
    CNST(227.),        CNST(232.0381),   CNST(231.03588),  CNST(238.0289),   &
    CNST(237.0482),    CNST(244.),       CNST(243.),       CNST(247.),       &
    CNST(247.),        CNST(251.),       CNST(254.),       CNST(257.),       &
    CNST(258.),        CNST(259.),       CNST(260.)                          &
   /)

  character(len=22), public :: pt_electronic_configuration(pt_n_elements) = (/ &
  "1s1                   ", "1s2                   ",  "[He] 2s1              ", "[He] 2s2              ",  &
  "[He] 2s2 2p1          ", "[He] 2s2 2p2          ",  "[He] 2s2 2p3          ", "[He] 2s2 2p4          ",  &
  "[He] 2s2 2p5          ", "[He] 2s2 2p6          ",  "[Ne] 3s1              ", "[Ne] 3s2              ",  &
  "[Ne] 3s1 3p1          ", "[Ne] 3s1 3p2          ",  "[Ne] 3s1 3p3          ", "[Ne] 3s1 3p4          ",  &
  "[Ne] 3s1 3p5          ", "[Ne] 3s1 3p6          ",  "[Ar] 4s1              ", "[Ar] 4s2              ",  &
  "[Ar] 3d1 4s2          ", "[Ar] 3d2 4s2          ",  "[Ar] 3d3 4s2          ", "[Ar] 3d5 4s1          ",  &
  "[Ar] 3d5 4s2          ", "[Ar] 3d6 4s2          ",  "[Ar] 3d7 4s2          ", "[Ar] 3d8 4s2          ",  &
  "[Ar] 3d10 4s1         ", "[Ar] 3d10 4s2         ",  "[Ar] 3d10 4s2 4p1     ", "[Ar] 3d10 4s2 4p2     ",  &
  "[Ar] 3d10 4s2 4p3     ", "[Ar] 3d10 4s2 4p4     ",  "[Ar] 3d10 4s2 4p5     ", "[Ar] 3d10 4s2 4p6     ",  &
  "[Kr] 5s1              ", "[Kr] 5s2              ",  "[Kr] 4d1 5s2          ", "[Kr] 4d2 5s2          ",  &
  "[Kr] 4d4 5s1          ", "[Kr] 4d5 5s1          ",  "[Kr] 4d5 5s2          ", "[Kr] 4d7 5s1          ",  &
  "[Kr] 4d8 5s1          ", "[Kr] 4d10             ",  "[Kr] 4d10 5s1         ", "[Kr] 4d10 5s2         ",  &
  "[Kr] 4d10 5s2 5p1     ", "[Kr] 4d10 5s2 5p2     ",  "[Kr] 4d10 5s2 5p3     ", "[Kr] 4d10 5s2 5p4     ",  &
  "[Kr] 4d10 5s2 5p5     ", "[Kr] 4d10 5s2 5p6     ",  "[Xe] 6s1              ", "[Xe] 6s2              ",  &
  "[Xe] 5d1 6s2          ", "[Xe] 4f2 6s2          ",  "[Xe] 4f3 6s2          ", "[Xe] 4f4 6s2          ",  &
  "[Xe] 4f5 6s2          ", "[Xe] 4f6 6s2          ",  "[Xe] 4f7 6s2          ", "[Xe] 4f7 5d1 6s2      ",  &
  "[Xe] 4f9 6s2          ", "[Xe] 4f10 6s2         ",  "[Xe] 4f11 6s2         ", "[Xe] 4f12 6s2         ",  &
  "[Xe] 4f13 6s2         ", "[Xe] 4f14 6s2         ",  "[Xe] 4f14 5d1 6s2     ", "[Xe] 4f14 5d2 6s2     ",  &
  "[Xe] 4f14 5d3 6s2     ", "[Xe] 4f14 5d4 6s2     ",  "[Xe] 4f14 5d5 6s2     ", "[Xe] 4f14 5d6 6s2     ",  &
  "[Xe] 4f14 5d7 6s2     ", "[Xe] 4f14 5d9 6s1     ",  "[Xe] 4f14 5d10 6s1    ", "[Xe] 4f14 5d10 6s2    ",  &
  "[Xe] 4f14 5d10 6s2 6p1", "[Xe] 4f14 5d10 6s2 6p2",  "[Xe] 4f14 5d10 6s2 6p3", "[Xe] 4f14 5d10 6s2 6p4",  &
  "[Xe] 4f14 5d10 6s2 6p5", "[Xe] 4f14 5d10 6s2 6p6",  "[Rn] 7s1              ", "[Rn] 7s2              ",  &
  "[Rn] 6d1 7s2          ", "[Rn] 6d2 7s2          ",  "[Rn] 5f2 6d1 7s2      ", "[Rn] 5f3 6d1 7s2      ",  &
  "[Rn] 5f4 6d1 7s2      ", "[Rn] 5f6 7s2          ",  "[Rn] 5f7 7s2          ", "[Rn] 5f7 6d1 7s2      ",  &
  "[Rn] 5f9 7s2          ", "[Rn] 5f10 7s2         ",  "[Rn] 5f11 7s2         ", "[Rn] 5f12 7s2         ",  &
  "[Rn] 5f13 7s2         ", "[Rn] 5f14 7s2         ",  "[Rn] 5f14 6d1 7s2     "                             &
  /)

contains
  integer function pt_number_from_symbol(s) result(z)
    character(len=3), intent(in)    :: s

    integer :: i

    z = -1
    do i = 1, pt_n_elements
      if(s == pt_symbol(i)) then
        z = i
        exit
      end if
    end do

  end function pt_number_from_symbol

end module periodic_table_m

!! Local Variables:
!! mode: f90
!! coding: utf-8
!! End:
