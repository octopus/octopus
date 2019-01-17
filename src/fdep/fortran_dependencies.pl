#!/usr/bin/perl -w
#
# Copyright 2015 Lorenz Hüdepohl
#
# This file is part of fdep and licensed under the MIT license
# see the file LICENSE for more information
#

use strict;

my %defs = ();
my %uses = ();
my %incs = ();
my %files = ();

# mode is the first argument: either mod or inc
my $mode = shift;

my $use_re = qr/^\s*use\s+(\S+)\s*$/;
my $def_re = qr/^\s*module\s+(\S+)\s*$/;
my $inc_re = qr/^\s*(\S+)\s*$/;

sub add_use {
	my ($file, $module) = @_;
	if (defined($defs{$module}) && $defs{$module} eq $file) {
		# do not add self-dependencies
		return;
	}
	if (!defined($uses{$file})) {
		$uses{$file} = { $module => 1 };
	} else {
		$uses{$file}{$module} = 1;
	}
}

sub add_def {
	my ($file, $module) = @_;
	if (!defined($defs{$module})) {
		$defs{$module} = $file;
		if (defined($uses{$file}) && defined($uses{$file}{$module})) {
			delete $uses{$file}{$module};
		}
	} else {
		die "Module $module both defined in $file, $defs{$module}";
	}
}

sub add_inc {
	my ($file, $module) = @_;
	if (!defined($incs{$file})) {
		$incs{$file} = { $module => 1 };
	} else {
		$incs{$file}{$module} = 1;
	}
}

my $target = shift;

foreach my $file (<>) {
	chomp($file);
	if (exists $files{$file}) {
		next;
	} else {
		$files{$file} = 1;
	}
	my $re;
	my $add;
	my $object;
        my $type;
	if (defined($ENV{V}) && $ENV{V} ge "2") {
		print STDERR "fdep: Considering file $file for target $target\n";
	}
	if ($file =~ /^.*\.\/.fortran_dependencies\/(.*)\.def_mods[^.]*(\..*)$/) {
		$re = $def_re;
		$add = \&add_def;
		$object = $1 . $2;
                $type = 'DEF'
	} elsif ($file =~ /^.*\.\/.fortran_dependencies\/(.*)\.use_mods[^.]*(\..*)$/) {
		$re = $use_re;
		$add = \&add_use;
		$object = $1 . $2;
                $type = 'USE'
	} elsif ($file =~ /^.*\.\/.fortran_dependencies\/(.*)\.inc_mods[^.]*(\..*)$/) {
		$re = $inc_re;
		$add = \&add_inc;
		$object = $1 . $2;
                $type = 'INC'
	} else {
		die "Unrecognized file extension for '$file'";
	}
	open(FILE,"<",$file) || die "\nCan't open $file: $!\n\n";
	while(<FILE>) {
		chomp;
                if ($type eq "DEF" or $type eq "USE") {
		        $_ = lc($_);
                }
		if ($_ =~ $re) {
			&$add($object, $1);
		} else {
			die "At $file:$.\nCannot parse module statement '$_', was expecting $re";
		}
	}
	close(FILE)
}

# module dependencies
if (lc($mode) eq 'mod') {
       foreach my $object (sort keys %uses) {
               for my $m (keys %{$uses{$object}}) {
                       if (defined $defs{$m}) {
                               print "$object: ", $defs{$m}, "\n";
                       } elsif (defined($ENV{V}) && $ENV{V} ge "1") {
                               print STDERR "Warning: Cannot find definition of module $m in files for current target $target, might be external\n";
                       }
               }
       }
# include file dependencies
} elsif (lc($mode) eq 'inc') {
        foreach my $object (sort keys %incs) {
                for my $m (keys %{$incs{$object}}) {
                        # This could be done in a more generic way to find the include files
                        print "$object: ", glob($m), "\n";
                }
        }
}
