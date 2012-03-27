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
# $Id$

use Getopt::Std;

getopts "hs:b:";

if($opt_h) {
    print <<"EndOfUsage";

Usage: var2html.pl [-b DIR] [-h]

    -b    The top level build tree directory, . if omitted
    -h    This help message

EndOfUsage

    exit 0;
}

$top_srcdir = ($opt_s ? $opt_s : ".");
$top_builddir = ($opt_b ? $opt_b : ".");

$doc = "$top_builddir/doc";
$share = "$top_builddir/share";

if(!-f "$share/varinfo_orig") {
    print stderr <<"EndOfErrorMsg";

The src and share directory could not be found. Please run
this script from the octopus toplevel directory or set -s and
-b options appropriately.

EndOfErrorMsg

    exit 1;
}

# configuration
$out_dir = "$doc/html/vars";


read_varinfo();   # generates %vars
sort_variables(); # generates @sorted

# make html pages
print_index();
print_alpha();
print_vars();
print_texi();

sub read_varinfo(){
  my $thisvar, $thisfield, $l;
  my $key, $arg;

  # let us parse the varinfo file
  open(IN, "<$share/varinfo_orig");
  $thisvar = "";
  $thisfield = "";
  while($l = <IN>){
    if($thisvar && $thisfield && $l !~ /^\w/){
      if($l =~ /^\s*$/){
	$vars{$thisvar}{$thisfield} .= "<br><br>\n";
      }else{
	$vars{$thisvar}{$thisfield} .= $l if($thisfield);
      }
      next;
    }else{
      $thisfield = "";
    }

    $l =~ /^(\S+)\s*(.*)$/;
    $key = lc($1);
    $arg = $2;

    if($key eq "variable"){
      $thisvar = $arg;
      $vars{$thisvar}{"variable"} = $thisvar;
    }

    next if(!$thisvar);

    if($key eq "section"){
      $vars{$thisvar}{"section"} = $arg;
    }elsif($key eq "type"){
      $vars{$thisvar}{"type"} = $arg;
    }elsif($key eq "default"){
      $vars{$thisvar}{"default"} = $arg;
    }elsif($key eq "description"){
      $thisfield = "description";
      $vars{$thisvar}{$thisfield} = "<br>";
    }elsif($key eq "option"){
      $arg =~ /^(\S*)\s*(.*?)\s*$/;
      $thisfield = "option_$2_$1";
    }elsif($key eq "end"){
      $thisvar   = "";
      $thisfield = "";
    }
  }
  close(IN);
}


sub sort_variables(){
  my %groups, @groups_sorted, $i;

  # we first make the groups
  foreach $key (keys %vars){
    $groups{$vars{$key}{"section"}} = 1;
  }
  @groups_sorted = sort keys %groups;

  for($i=0; $i<=$#groups_sorted; $i++){
    my @new = ();
    foreach $key (keys %vars){
      push @new, $key if($vars{$key}{"section"} eq $groups_sorted[$i]);
    }
    push @sorted, sort @new;
  }
}


sub print_index(){
  my $i, $j, $k, $key, $sect, $new_sect, $level;
  my @old, @new;

  open(OUT_index, ">$out_dir/vars_index.html");
  open(OUT_tree,  ">$out_dir/sections.js");

  print OUT_index "<ul>\n";
  print OUT_tree  "
USETEXTLINKS = 1
STARTALLOPEN = 0
ICONPATH = 'icons/'

aux0 = gFld(\"Variables\", \"varsRightFrame.html\")
foldersTree = aux0
";

  $sect = "";
  $level = 0;
  for($i=0; $i<=$#sorted; $i++){
    $key = $sorted[$i];

    $new_sect = $vars{$key}{"section"};
    if($sect ne $new_sect){
      @old = split("::", $sect);
      @new = split("::", $new_sect);

      for($j=0; $old[$j] && $new[$j] && $old[$j] eq $new[$j]; $j++){};

      for($k=$j; $old[$k]; $k++){
	print OUT_index "</ul>\n";
      }

      for($k=$j; $new[$k]; $k++){
	$new_sect =~ /^([^:]*)/;
	$link     = "$1.html#".$vars{$key}{"section"};
	$link     =~ s/\s/_/g;

	print OUT_index "<li>", $new[$k], "</li>\n<ul>\n";
	print OUT_tree  "aux", $k+1, " = insFld(aux$k, gFld(\"", $new[$k], "\", \"vars/$link\"))\n";
      }

      $sect = $new_sect;
      $level = $k;
    }

    $sect =~ /^([^:]*)/;
    $link = "$1.html#".$vars{$key}{"variable"};
    $link =~ s/\s/_/g;

    print OUT_index "<li><a href='$link'>", $vars{$key}{"variable"}, " </a></li>\n";
    print OUT_tree  "insDoc(aux$level, gLnk(\"R\", \"", $vars{$key}{"variable"}, "\", \"vars/$link\"))\n";
  }

  #print last <\ul>
  @old = split("::", $sect);
  for($j=0; $old[$j]; $j++){
    print OUT_index "</ul>\n";
  }
  print OUT_index "</ul>\n";

  # close files
  close(OUT_index);
  close(OUT_tree);
}


sub print_alpha(){
  my $key, $letter;

  open(OUT_tree,  ">$out_dir/alpha.js");

  print OUT_index "<ul>\n";
  print OUT_tree  "
USETEXTLINKS = 1
STARTALLOPEN = 0
ICONPATH = 'icons/'

aux0 = gFld(\"Variables\", \"varsRightFrame.html\")
foldersTree = aux0
";

  $letter = '';
  foreach $key (sort { lc($a) cmp lc($b) } keys %vars){
    $key =~ /^(.)/;
    $new_letter = lc($1);

    if($new_letter ne $letter){
      $letter = $new_letter;
      print OUT_tree  "aux1 = insFld(aux0, gFld(\"", uc($letter), "\", \"javascript:parent.op()\"))\n";
    }

    $vars{$key}{"section"} =~ /^([^:]*)/;
    my $link = "vars/$1.html#".$vars{$key}{"variable"};
    $link =~ s/\s/_/g;
    print OUT_tree  "insDoc(aux1, gLnk(\"R\", \"", $vars{$key}{"variable"}, "\", \"$link\"))\n"; 
  }

  close(OUT_tree);
}


sub print_vars(){
  my $i, $k, $sect, $new_sect;

  $sect = "";
  $sect_full = "";
  for($i=0; $i<=$#sorted; $i++){
    $key = $sorted[$i];

    $new_sect = $vars{$key}{"section"};
    $new_sect =~ s/:.*$//;
    $new_sect =~ s/\s/_/g;

    if($sect ne $new_sect){
      my $link;

      if($sect){
	print OUT "
</body>
</html>";
	close(OUT)
      }
      $sect = $new_sect;

      open(OUT, ">$out_dir/$sect.html");
      print OUT "
<html>
<head>

<style>
   BODY {background-color: white; 
         font-size: 10pt; font-family: verdana,helvetica;}
   A  {text-decoration: none;color: blue}
</style>
</head>
<body>
";
    }

    if($sect_full ne $vars{$key}{"section"}){
      $sect_full = $vars{$key}{"section"};
      print OUT "\n<a name='", $vars{$key}{"section"}, "'</a>\n";
      print OUT "<H2>", $vars{$key}{"section"}, "</H2>\n";
    }

    # print the information about the variable
    print OUT
      "\n
<p><b><a name='$vars{$key}{variable}'></a>$vars{$key}{variable}</b>
<br/><i>Section</i>: $vars{$key}{section}
<br/><i>Type</i>: $vars{$key}{type}\n";
    print OUT "<br/><i>Default</i>: $vars{$key}{default}\n" if($vars{$key}{default});
    print OUT "<br/>$vars{$key}{description}\n";

    my $first = 1;
    foreach $k (sort keys %{ $vars{$key} }) {
      if($k =~ /^option_(.*?)_(.*)/) {
	if($first){
	  print OUT "<br/><i>Options</i>:\n<ul>\n";
	  $first = 0;
	}
	print OUT "<li><b>$2</b>";
	print OUT ": ", $vars{$key}{$k}, "</li>\n";
      }
    }
    print OUT "</ul>\n" if(!$first);
    print OUT "</p><hr width='30%' align='left'/>\n";
  }
  close(OUT);
}


sub print_texi(){
  my $i, $k, $sect, $new_sect;

  open(OUT, ">$doc/variables.texi");
  print OUT "\@node Input Variables,,,
\@chapter Input Variables

\@code{octopus} has quite a few options, that we will subdivide in different groups.
After the name of the option, its type and default value (when applicable)
are given in parenthesis.
";

  $sect = "";
  $sect_full = "";
  for($i=0; $i<=$#sorted; $i++){
    $key = $sorted[$i];

    $new_sect = $vars{$key}{"section"};
    $new_sect =~ s/:.*$//;
    $new_sect =~ s/\s/_/g;

    if($sect_full ne $vars{$key}{"section"}){
      my $i, @p;

      if($sect_full){
	print OUT "\@end itemize\n";
      }

      $sect_full = $vars{$key}{"section"};

      @p = split("::", $sect_full);
      print OUT '@node ', to_texi($p[$#p]), ",,,\n";
      print OUT '@';
      for($i=0; $i<$#p; $i++){
	print OUT "sub";
      }
      print OUT "section ", to_texi($p[$#p]), "\n";
      print OUT "\@c ----------------------------------\n\n";
      print OUT "\@itemize\n";
    }

    # print the information about the variable
    print OUT "\@item \@strong{", to_texi($vars{$key}{variable}), "}@*\n",
      "\@vindex \@code{", to_texi($vars{$key}{variable}), "}@*\n",
      "\@emph{Section}: ", to_texi($vars{$key}{section}), "@*\n",
      "\@emph{Type}: ", to_texi($vars{$key}{type}), "@*\n";
    print OUT "\@emph{Default}: ", to_texi($vars{$key}{default}), "@*\n" if($vars{$key}{default});
    print OUT to_texi($vars{$key}{description}), "\n\n";

    my $first = 1;
    foreach $k (sort keys %{ $vars{$key} }) {
      if($k =~ /^option_(.*?)_(.*)/) {
	if($first){
	  print OUT "\@emph{Options}:\n\@itemize \@minus\n";
	  $first = 0;
	}
	print OUT "\@item \@strong{", to_texi($2), "}";
	print OUT ": ", to_texi($vars{$key}{$k});
      }
    }
    print OUT "\@end itemize\n" if(!$first);

    print OUT "\n\@c ----------------------------------\n";
  }

  # close the last itemize
  print OUT "\@end itemize\n";

  # close file
  close(OUT);
}

sub to_texi(){
  my ($t) = @_;

  #$t =~ s#\{#\@\{#g;
  #$t =~ s#\}#\@\}#g;

  $t =~ s#<br/*><br/*>##gi;
  $t =~ s#<br/*>#@*#gi;
  $t =~ s#<hr/*>#------------------------------------------@*#gi;
  $t =~ s#&nbsp;#@ #gi;

  $t =~ s#<i>#\@emph\{#gi;
  $t =~ s#<b>#\@strong\{#gi;
  $t =~ s#</[bi]>#\}#gi;

  while($t =~ m#<MATH>(.*?)</MATH>#s){
    $tex      = $1;
    $tex      =~ s#\n# #gs;

    $verbatim = $1;
    $verbatim =~ s#\{#\@\{#g;
    $verbatim =~ s#\}#\@\}#g;
    $verbatim =~ s#\n# #gs;

    $t =~ s#<MATH>(.*?)</MATH>#\@ifnottex\n\@verbatim\n$verbatim\n\@end verbatim\n\@end ifnottex\n\@tex\n\$\$\n$tex\n\$\$\n\@end tex\n#s;
  }
  $t =~ s#<math>(.*?)</math>#\@math\{$1\}#gs;

  $t =~ s#<tt>#\@t\{#gi;
  $t =~ s#</tt>#\}#gi;

  $t =~ s#<ul>#\@itemize\n#gi;
  $t =~ s#</ul>#\@end itemize\n#gi;
  $t =~ s#<li>#\@item\n#gi;
  $t =~ s#</li>##gi;

  $t =~ s#&gt;#>#gi;
  $t =~ s#&lt;#<#gi;

  return $t;
}
