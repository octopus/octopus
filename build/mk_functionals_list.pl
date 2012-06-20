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
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
# 02111-1307, USA.
#
# $Id: mk_varinfo.pl 5657 2009-06-29 15:02:07Z marques $


use Getopt::Std;
use File::Find;
getopts "hs:b:";

if($opt_h) {
    print <<"EndOfUsage";

Usage: mk_functionals_list.pl [-b DIR] [-s DIR] [-h]

    -b    The top level build tree directory, . if omitted
    -s    The top level source tree directory, . if omitted
    -h    This help message
EndOfUsage

    exit 0;
}

$top_srcdir = ($opt_s ? $opt_s : ".");
$top_builddir = ($opt_b ? $opt_b : ".");

$src   = "$top_srcdir/src/xc";
$funct = "$top_builddir/libxc/src/xc_funcs.h";

if(!-d $src && !-f $funct) {
    print stderr <<"EndOfErrorMsg";

The $src directory or the file $funct could not be found. Please run
this script from the octopus toplevel directory or set -s option appropriately.
EndOfErrorMsg
}

open(OUT, ">$src/functionals_list.F90");
print OUT <<"EndOfHeader";
! Note: this file is generated automatically by build/mk_functionals_list.pl
!
!%Variable XCFunctional
!%Type integer
!%Section Hamiltonian::XC
!%Description
!% Defines the exchange and correlation functional to be used;
!% they should be specified as a sum of a correlation term and an
!% exchange term. Defaults:
!% <br> 1D: lda_x_1d + lda_c_1d_csc
!% <br> 2D: lda_x_2d + lda_c_2d_amgb
!% <br> 3D: lda_x + lda_c_pz_mod
EndOfHeader

open(IN, "<$funct");
while($_ = <IN>){
  if(/\#define\s+(\S+)\s+(\d+)\s*\/\*\s*(.*?)\s*\*\//){
    $option  = $1;
    $number  = $2;
    $comment = $3;

    if($option =~ /^XC_\S+_C_/ || $option =~ /^XC_\S+_XC_/){
      $number *= 1000;
    }

    $option =~ s/XC_(.*)$/\L$1\E/g;
    print OUT "!%Option $option               $number\n!% $comment\n";
  }
}
print OUT <<EOF;
!%Option oep_x                    901
!% OEP: Exact exchange
!%Option lda_xc_cmplx             701
!% LDA complex scaled exchange-correlation.
!%Option none                       0
!% Exchange and correlation set to zero.
!%End
EOF

close(IN);
close(OUT);

