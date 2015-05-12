#!/usr/bin/env perl

# Copyright (C) 2013-2015 D. Strubbe
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
# $Id$

# Query a buildbot using the Octopus/APE/BerkeleyGW testsuite infrastructure for
# values obtained by each buildslave for a particular test and match.
# Parses the HTML status pages.
# Tested with BuildBot 0.7.12, 0.8.5, 0.8.8

if(@ARGV < 2) {
  print "Usage: buildbot_query.pl <inputfile> <match> [<branch> [localfiles]]\n";
  exit(1);
}

$inputfile = shift @ARGV;
$match = shift @ARGV;
$branch = shift @ARGV;
if(length($branch) == 0) {
    $branch = "trunk";
}

# options specifying setup for Octopus
$bbpath = "http://www.tddft.org/programs/octopus/buildbot";
$shell_num = 3;

print "URL: $bbpath\n";
print "Branch: $branch\n";
print "Input file: $inputfile\n";
print "Match: $match\n\n";

$branch = "?branch=$branch";

if($inputfile =~ ".test") {
    print STDERR "Pass the input file, not the .test file.\n";
    exit(1);
}

# get list of latest builds
if(-e "one_box_per_builder") { system ("rm one_box_per_builder"); }
system("wget -nv $bbpath/builders$branch -O one_box_per_builder");

if(@ARGV > 0) {
    open(ONEBOX, ">>one_box_per_builder") or die "Cannot open one_box_per_builder.\n";
    print ONEBOX "\n";
    foreach(@ARGV) {
	print ONEBOX "LOCAL FILENAME: $_\n";
    }
    close(ONEBOX);
}

open(ONEBOX, "<one_box_per_builder") or die "Cannot open one_box_per_builder.\n";

$total = 0.0;
$counts = 0;
$max = -inf;
$maxname = "";
$min = inf;
$minname = "";

$filename = "text";

while ($_ = <ONEBOX>) {
# BB 0.7.12
# <td align="center" class="LastBuild box success"><a href="builders/lascar_x86_64_gfortran_cl_intel/builds/139">10187</a><br />build<br />successful</td>
# BB 0.8.5
#<a href="builders/mauchly_x86_64_intel_openmp/builds/80">10898</a>
# BB 0.8.8
#<td class="box"><a href="./builders/mauchly_x86_64_intel_openmp">mauchly_x86_64_intel_openmp</a></td>
    if ( $_ =~ /<a href=".*builders\/(.*)\/builds\/(.*)">(.*)<\/a>/) {
	$builder = $1;
	$build_num = $2;
	$svn_rev = $3;
	# rebuild the URL
	$url = "builders/$builder/builds/$build_num";
	print "\nBuilder: $builder, at svn revision $svn_rev\n";

	# remove old file, or new ones will be named 'text.2' etc.
	if(-e $filename) { system ("rm $filename"); }
	# get the version without any complicating HTML markup
	unless (system ("wget -nv $bbpath/$url/steps/shell_$shell_num/logs/stdio/$filename") == 0) {
	    print STDERR "Cannot download test log.\n";
	    next;
	}

	$name = $builder;
	unless (open(TESTLOG, "<$filename")) {
	    print STDERR "Cannot open test log.\n";
	    next;
	}
    } elsif ( $_ =~ /LOCAL FILENAME: (.*)/) {
	print "\n$_\n";
	$name = $1;
	open(TESTLOG, "<$1");
    } else {
	next;
    }
    $file_found = 0;
    $match_found = 0;
    while ($_ = <TESTLOG>) {
	
	# do not use ~= / .. / here or $inputfile needs to have special characters escaped
	if($_ =~ /Using input file : / && index($_, $inputfile) != -1) {
	    $file_found = 1;
	    while ($_ = <TESTLOG>) {
		if(index($_, $match) != -1) {
		    if($_ =~ /\(Calculated value = (.*)\)/) {  # match OK
			print $_;
			$value = $1;
			$match_found = 1;
		    } else {  # match FAIL
			while ($_ = <TESTLOG>) {
			    if(index($_, $match) == -1) {
				if($_ !~ /^$/) { # print if not blank
				    print $_;
				    if($_ =~ /Calculated value : (.*)/) {
					$value = $1;
					$match_found = 1;
				    }
				}
			    } else {
				last;
			    }
			}
		    }
		    
		    # If match failed and did not give a number, do not treat it as zero.
		    if($match_found) {
			$total += $value;
			$counts += 1;
			if($value < $min) {
			    $minname = $name;
			    $min = $value;
			}
			if($value > $max) {
			    $maxname = $name;
			    $max = $value;
			}
		    }
		}
		if($_ =~ /Using input file/) { last; }
	    }
	}
    }
    close(TESTLOG);
    if($file_found == 0) {
	print STDERR "File not found.\n";
    } elsif($match_found == 0) {
	print STDERR "Match not found.\n";
    }
    # why not? builder down, svn or compilation failed, not in the right category of builders, etc.
}

if($counts == 0) {
    print STDERR "No matches found.\n";
    exit(1);
} else {
    print "\n\n=== SUMMARY ===\n";
    print "Based on $counts matches found.\n";
    print "Minimum   = $min\n";
    print "    ($minname)\n";
    print "Maximum   = $max\n";
    print "    ($maxname)\n";
    print "Average   = " . ($total / $counts) . "\n\n";
    print "Center    = " . ($max + $min)/2 . "\n";
    printf "Precision = %e\n", ($max - $min)/2;
}

close(ONEBOX);
