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
getopts "hs:b:";

if($opt_h) {
    print <<"EndOfUsage";

Usage: mk_varinfo.pl [-b DIR] [-s DIR] [-h]

    -b    The top level build tree directory, . if omitted
    -s    The top level source tree directory, . if omited
    -h    This help message

EndOfUsage

    exit 0;
}

$top_srcdir = ($opt_s ? $opt_s : ".");
$top_builddir = ($opt_b ? $opt_b : ".");

$src = "$top_srcdir/src";
$share = "$top_builddir/share";
$include = "$top_builddir/src/include";

if(!-d $src && !-d $share) {
    print stderr <<"EndOfErrorMsg";

The src and share directory could not be found. Please run
this script from the octopus toplevel directory or set -s and
-b options appropriately.

EndOfErrorMsg

    exit 1;
}

# get all files named *.F90 recursively
@F90 = ();
finddepth (sub{ push @F90, $File::Find::name if /\.F90$/ }, "$src");

# Abort with warning if no *.F90 files were found.
if($#F90 < 0) {
    print stderr <<"EndOfWarning";

Warning: No *.F90 files found. Probably, the source directory
was not set correctly.

EndOfWarning

   exit 2;
}

open(OUT_text, ">$share/varinfo");
open(OUT_orig, ">$share/varinfo_orig");

%opt_value = ();
%varopt = ();
%default_value = ();

foreach $F90file (@F90){
  open(IN, "<$F90file");
  while($_=<IN>){
    if(/!%Variable\s+(\S+)/i){
      $var = $1;

      $desc = "";
      do {
	s/^\s*!%//; s/\s*$//;

	if(/^Default\s+(\S+)/){
	  put_default($1, $var);
	}
	    
	if(/^Option\s+(\S+)\s+bit\((\S+)\)/){
	  if($2 > 52) {
	    printf STDERR "ERROR: bit($2) is too large and will overflow the maximum integer.\n";
	    printf STDERR "File $F90file, Variable $var, Option $1.\n";
	    exit(1);
	  }
	  put_opt($1, (1<<($2)), $var);
	} elsif(/^Option\s+(\S+)\s+(\S+)/){
	  put_opt($1, $2, $var);
	}

	if(/^ </){ # lines that start with a html command
	  print_desc($desc) if($desc ne "");
	  $desc = $_."\n";
	}elsif(/^ /){ # get whole description
	  $desc .= $_."\n";
	}else{
	  if($desc ne "") {
	    print_desc($desc);
	    $desc = "";
	  }

	  printf OUT_text "%s\n", $_;
	  printf OUT_orig "%s\n", $_;
	}

	$_ = <IN>;
      } until (eof(IN) || !/!%/ || /!%End/i);

      if($desc ne "") {
	print_desc($desc);
      }

      if(!/!%End/i){
	print stderr "In info block of variable $var (file src/$F90file), block is incomplete\n";
	exit 2;
      }else{
	printf OUT_text "END\n\n";
	printf OUT_orig "END\n\n";
      }
    }
  }
  close(IN);
}

close(OUT_text);
close(OUT_orig);

print_opt();
print_options_header();
print_defaults_header();

#####################################################
# tries to put an option in global %opt_value
sub put_opt{
  my ($a, $b, $var) = @_;

  if($opt_value{$a} && ($opt_value{$a} ne $b)){
    print stderr "Option '", $a, "' is multiply defined\n",
      "    '", $opt_value{$b}, "' ne '", $b, "'\n";
    exit 3;
  }
  
  $opt_value{$a} = $b;
  $varopt{$var."__".$a} = $b;
 
}

#####################################################
# tries to put an option in global %opt_value
sub put_default{
  my ($a, $var) = @_;

  $default_value{$var} = $a;
 
}

#####################################################
# prints %opt_value to share/variables
sub print_opt{
  open(OUT, ">$share/variables");
  my $key;
  foreach $key (sort(keys %opt_value)) {
    print OUT $key, " = ", $opt_value{"$key"}, "\n";
  }

  close(OUT);
}

#####################################################
# prints %default_value to src/include/defaults.h
sub print_defaults_header{
  open(OUT, ">$include/defaults.h");
  my $key;
  foreach $key (sort(keys %default_value)) {
    print OUT "#define DEFAULT__", uc $key, " (", $default_value{"$key"}, "_8)\n";
  }
  close(OUT);
}

#####################################################
# prints %opt_value to src/include/options.h
sub print_options_header{
  open(OUT, ">$include/options.h");
  my $key;
  foreach $key (sort(keys %varopt)) {
    print OUT "#define OPTION__", uc $key, " (", $varopt{"$key"}, "_8)\n";
  }

  close(OUT);
}


#####################################################
# justifies a string
sub print_desc(){
  my ($desc) = @_;
  my $ml = 75;

  my $i, $line;

  print OUT_orig $desc;

  $desc =~ s/\n/ /gm;
  $desc =~ s/^\s*//;
  $desc =~ s/\s\s+/ /g;

  # convert html to text
  $desc =~ s#<br/*>##g;
  $desc =~ s#<hr/*>#------------------------------------------\n#g;
  $desc =~ s#&nbsp;# #g;

  $desc =~ s#</*i>#_#g;
  $desc =~ s#</*b>#*#g;
  $desc =~ s#</*math>##g;
  $desc =~ s#</*tt>##g;

  $desc =~ s#</*ul>##g;
  $desc =~ s#<li># *) #g;
  $desc =~ s#</li>##g;

  while(){
    if(length($desc) <= $ml){
      print OUT_text " ", $desc, "\n";
      last;
    }

    for($i=$ml; $i>0 && substr($desc, $i, 1) ne ' '; $i--){};
    if($i == 0){
      for($i=0; $i<length($desc) && substr($desc, $i, 1) ne ' '; $i++){};
    }

    $line = substr($desc, 0, $i);
    $desc = substr($desc, $i+1);

    $spaces = $ml - length($line);
    if($spaces > 0){
      for($i=0; $i<$spaces; $i++){
	$line =~ s/([^ ]) ([^ ])/$1  $2/;
      }
    }
    print OUT_text " ", $line, "\n";
  }
}
