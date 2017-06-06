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

Usage: mk_kfunctionals_list.pl [-s DIR] [-I DIR] [-h]

    -s    The top-level source tree directory, . if omitted
    -I    The libxc include directory, in which to find xc_funcs.h
    -h    This help message
EndOfUsage

    exit 0;
}

$top_srcdir = ($opt_s ? $opt_s : ".");

$src   = "$top_srcdir/src/frozen";
$vrsn  = "$opt_I/xc_version.h";
$funct = "$opt_I/xc_funcs.h";


if(!-d $src) {
    print STDERR "Cannot find directory '$src'. Run from top-level directory or set -s option appropriately.\n";
    exit(1);
}

if(!-f $vrsn) {
    print STDERR "Cannot find file '$vrsn'. Set -I option appropriately.\n";
    exit(1);
}

if(!-f $funct) {
    print STDERR "Cannot find file '$funct'. Set -I option appropriately.\n";
    exit(1);
}

open(IN, "<$vrsn");
while($_ = <IN>){
    if(/\#define\s+XC\_VERSION\s+\"(\S+)\"/){
	$version  = $1;
	last;
    }
}
close(IN);

open(OUT, ">$src/kfunctionals_list.F90");
print OUT <<"EndOfHeader";
! Note: this file is generated automatically by build/mk_kfunctionals_list.pl
!
!%Variable TnaddFunctional
!%Type integer
!%Section Hamiltonian::Subsystems
!%Description
!% Defines the Kinetic Functional to be used in a Subsystem calculation,
!% For more information on the functionals, see
!% <a href=http://octopus-code.org/wiki/Libxc:manual#Available_functionals>
!% Libxc documentation</a>. The list provided here is from libxc $version; if you have
!% linked against a different libxc version, you may have a somewhat different set
!% of available functionals.
!% <br>Default: <tt>none</tt>
EndOfHeader

open(IN, "<$funct");
while($_ = <IN>){
    if(/\#define\s+(\S+)\s+(\d+)\s*\/\*\s*(.*?)\s*\*\//){
	$option  = $1;
	$number  = $2;
	$comment = $3;

	# include only kinetic-energy functionals
	next unless($option =~ /^XC_\S+_K_/);

	$option =~ s/XC_(.*)$/\L$1\E/g;
	print OUT "!%Option $option               $number\n!% $comment\n";
    }
}
print OUT <<EOF;
!%Option none                       0
!% Exchange and correlation set to zero (not from libxc).
!%End
EOF

close(IN);
close(OUT);
