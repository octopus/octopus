#!/usr/bin/env perl
#
# Copyright (C) 2002-2006 M. Marques, A. Castro, A. Rubio, G. Bertsch
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA.
#


use Getopt::Std;
use File::Find;
getopts "hs:I:";

if($opt_h) {
    print <<"EndOfUsage";

Usage: mk_functionals_list.pl [-s DIR] [-I DIR] [-h]

    -s    The top-level source tree directory, . if omitted
    -I    The libxc include directory, in which to find xc_funcs.h
    -h    This help message
EndOfUsage

    exit 0;
}

$top_srcdir = ($opt_s ? $opt_s : ".");

$src   = "$top_srcdir/src/hamiltonian";
$funct = "$opt_I/xc_funcs.h";

if(!-d $src) {
    print STDERR "Cannot find directory '$src'. Run from top-level directory or set -s option appropriately.\n";
    exit(1);
}

if(!-f $funct) {
    print STDERR "Cannot find file '$funct'. Set -I option appropriately.\n";
    exit(1);
}

open(OUT, ">$src/functionals_list.F90");
print OUT <<"EndOfHeader";
! Note: this file is generated automatically by build/mk_functionals_list.pl
!
!%Variable XCFunctional
!%Type integer
!%Section Hamiltonian::XC
!%Description
!% Defines the exchange and correlation functionals to be used,
!% specified as a sum of an exchange functional and a
!% correlation functional, or a single exchange-correlation functional
!% (<i>e.g.</i> <tt>hyb_gga_xc_pbeh</tt>). For more information on the functionals, see
!% <a href=https://gitlab.com/libxc/libxc/wikis/Functionals-list-4.0.4>
!% Libxc documentation</a>. The list provided here is from libxc 4; if you have
!% linked against a different libxc version, you may have a somewhat different set
!% of available functionals. Note that kinetic-energy functionals are not supported.
!%
!% The default functional will be selected by Octopus to be consistent
!% with the pseudopotentials you are using. If you are not using
!% pseudopotentials, Octopus cannot determine the functional used to
!% generate the pseudopotential, or the pseudopotential functionals
!% are inconsistent, Octopus will use the following defaults:
!%
!% <br>1D: <tt>lda_x_1d + lda_c_1d_csc</tt>
!% <br>2D: <tt>lda_x_2d + lda_c_2d_amgb</tt>
!% <br>3D: <tt>lda_x + lda_c_pz_mod</tt>
EndOfHeader

open(IN, "<$funct");
while($_ = <IN>){
  if(/\#define\s+(\S+)\s+(\d+)\s*\/\*\s*(.*?)\s*\*\//){
    $option  = $1;
    $number  = $2;
    $comment = $3;

    # do not include kinetic-energy functionals
    next if($option =~ /^XC_\S+_K_/);

    if($option =~ /^XC_\S+_C_/ || $option =~ /^XC_\S+_XC_/){
      $number *= 1000;
    }

    $option =~ s/XC_(.*)$/\L$1\E/g;
    print OUT "!%Option $option               $number\n!% $comment\n";
  }
}
print OUT <<EOF;
!%Option oep_x                    901
!% OEP: Exact exchange (not from libxc).
!%Option ks_inversion             801
!% Inversion of KS potential (not from libxc).
!%Option lda_xc_cmplx             701
!% Complex-scaled LDA exchange and correlation (not from libxc).
!%Option pbe_xc_cmplx             702
!% Complex-scaled PBE exchange and correlation (not from libxc).
!%Option lb94_xc_cmplx            703
!% Complex-scaled LB94 exchange and correlation (not from libxc).
!%Option rdmft_xc_m               601
!% RDMFT Mueller functional (not from libxc).
!%Option xc_half_hartree          917
!% Half-Hartree exchange for two electrons (supports complex scaling) (not from libxc).
!% Defined by <math>v_{xc}(r) = v_H(r) / 2</math>.
!%Option hyb_gga_xc_mvorb_hse06   921000 
!% Density-based mixing parameter of HSE06 (not from libxc).
!%Option hyb_gga_xc_mvorb_pbeh    922000
!% Density-based mixing parameter of PBEH (not from libxc).
!%Option vdw_c_vdwdf      918000
!% van der Waals density functional vdW-DF correlation from libvdwxc (not from libxc).  Use with gga_x_pbe_r.
!%Option vdw_c_vdwdf2     919000
!% van der Waals density functional vdW-DF2 correlation from libvdwxc (not from libxc).  Use with gga_x_rpw86.
!%Option vdw_c_vdwdfcx    920000
!% van der Waals density functional vdW-DF-cx correlation from libvdwxc (not from libxc).  Use with gga_x_lv_rpw86.
!%Option none                       0
!% Exchange and correlation set to zero (not from libxc).
!%End
EOF

close(IN);
close(OUT);
