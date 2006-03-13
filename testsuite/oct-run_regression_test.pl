#!/usr/bin/perl
#
# $Id$

use Getopt::Std;
use File::Basename;


sub usage {

  print <<EndOfUsage;

 Copyright (C) 2005 by Heiko Appel

Usage: oct-run_regression_test.pl [options]

    -n        dry-run
    -v        verbose
    -h        this usage
    -e        name of octopus executable
    -c        create template
    -f        filename of testsuite
    -d        working directory for the tests
    -i        print inputfile
    -p        preserve working directories
    -l        copy output log to current directory
    -m        run matches only (assumes there are work directories)

Report bugs to <appel\@physik.fu-berlin.de>.
EndOfUsage
  exit 0;
}


sub create_template {
  $date = `date +"%d.%m.%y"`;
  chomp($date);
  $arch = `uname -a`;
  chomp($arch);
  $author = `whoami`;
  chomp($author);
  $author =~ s/^(\w)(.*)/\u$1$2/;

  open(TEMPLATE, ">".$opt_c );

  print TEMPLATE <<EndOfTemplate;
# \$Id$

Test     :
Author   : $author
Date     : $date
Arch     : $arch
Release  :
Programs : octopus
Runtime  : short-run
Enabled  : Yes

Input: 01-template.01-ground_state.inp

RUN

# add your own matches
#
# Examples (of course you have to uncomment the lines :)
# Match ; TotalEnergy ; grep -A 6 '^Energy:' static/info ; \\s*Total\\s*=\\s*-0.33802679
# Match ; Eigenvalues ; grep -A 2 'Eigenvalues \\[' static/info | tail -1 ; \\s*-0.13735\\d*\\s*2.0000\\d*
EndOfTemplate

  close(TEMPLATE);

  print "Template written to: $opt_c \n";

  exit 0;

}

if (not @ARGV) { usage; }

getopts("nlvhe:c:f:d:ipm");

# Default values
$workdir = `mktemp -d /tmp/octopus.XXXXXX`;
chomp($workdir);

# Handle options
$opt_h && usage;
$opt_c && create_template;
if($opt_d)  { $workdir = $opt_d; }

if($opt_e)  { $octopus_exe = $opt_e; }
else{
  $octopus_exe =        "octopus";
  $octopus_exe = "../src/octopus"             if( -x "../src/octopus");
  $octopus_exe = "../src/octopus_mpi"         if( -x "../src/octopus_mpi");
  $octopus_exe = "../src/octopus_cmplx"       if( -x "../src/octopus_cmplx");
  $octopus_exe = "../src/octopus_cmplx_mpi"   if( -x "../src/octopus_cmplx_mpi");
}
if($octopus_exe !~ /^\//){
  $octopus_exe = $ENV{PWD}."/$octopus_exe";
}

die "could not find executable: $octopus_exe \n" if (`which $octopus_exe` eq "");


chomp($octopus_exe);

$octopus_base = basename($octopus_exe);
chomp($octopus_base);

# MPI stuff
$mpirun = $ENV{MPIRUN};
if ("$mpirun" eq "") { $mpirun = `which mpirun`; }
chomp($mpirun);
# default number of processors for MPI runs is 2
$np = 2;


if (!$opt_m) {
  system ("rm -rf $workdir");
  mkdir $workdir;
}

open(TESTSUITE, "<".$opt_f ) or die "cannot open testsuite file \n";


while ($_ = <TESTSUITE>) {

 # skip comments
 next if /^#/;

 if ( $_ =~ /Test\s*:\s*(.*)\s*$/) {
  $test{"name"} = $1;
  if(!$opt_i) {
    print "\033[34m ***** $test{\"name\"} ***** \033[0m \n\n";
    print "Using workdir    : $workdir \n";
    print "Using executable : $octopus_exe \n";
    print "Using test file  : $opt_f \n";
  }
 }


 if ( $_ =~ /Programs\s*:\s*(.*)\s*$/) {
  $valid=1;
  foreach my $program (split(/;/,$1)) {
    $program =~ s/^\s*//;
    $program =~ s/\s*$//;
    $valid = $program cmp $octopus_base;
    last if ( ! $valid );
    $program = "${program}_debug";
    $valid = $program cmp $octopus_base;
    last if ( ! $valid );
  }
  if ( $valid ) {
    print "\n$opt_f \033[31mnot valid\033[0m for executable: $octopus_exe\n";
    print "Skipping ... \n\n";
    exit 1;
  }
 }

 if ( $_ =~ /^Enabled\s*:\s*(.*)\s*$/) {
  %test = ();
  $enabled = $1;
  $enabled =~ s/^\s*//;
  $enabled =~ s/\s*$//;
  $test{"enabled"} = $enabled;
 }


 # Running this regression test if it is enabled
 if ( $enabled eq "Yes" ) {

   if ( $_ =~ /^Processors\s*:\s*(.*)\s*$/) {
     $np = $1;
   }

   if ( $_ =~ /^Input\s*:\s*(.*)\s*$/) {
     $input_file = dirname($opt_f) . "/" . $1;
     if( -f $input_file ) {
       print "\n\nUsing input file : $input_file \n";
       system("cp $input_file $workdir/inp");
     } else {
       die "could not find input file: $input_file\n";
     }

     if ( !$opt_m ) {
       if ( !$opt_n ) {
	 print "\nStarting test run ...\n";

	 # serial or MPI run?
	 if ( $octopus_exe =~ /mpi$/) {
	   if( -x "$mpirun") {
	     print "Executing: cd $workdir; $mpirun -np $np $octopus_exe &> out \n";
	     system("cd $workdir; $mpirun -np $np $octopus_exe &> out");
	   } else {
	     print "No mpirun found: Skipping parallel test \n";
	     exit 1;
	   }
	 } else {
	   print "Executing: cd $workdir; $octopus_exe &> out \n";
	   system("cd $workdir; $octopus_exe &> out");
	 }
	 system("grep -B2 -A5 'Running octopus' $workdir/out > build-stamp");
	 print "Finished test run.\n\n"; }
       else {
	 if(!$opt_i) { print "cd $workdir; $octopus_exe < inp &> out \n"; }
       }
       $test{"run"} = 1;
     }

   }

   if ( $_ =~ /^Match/ && !$opt_n) {
     $- = s/\\;/_COLUMN_/g;
     ($match, $name, $pre_command, $regexp) = split(/;/, $_);
     $name        =~ s/_COLUMN_/;/g;
     $pre_command =~ s/_COLUMN_/;/g;
     $regexp      =~ s/_COLUMN_/;/g;

     if(!opt_m) {
       die "have to run before matching" if !$test{"run"};
     }

     if ($opt_v) { print "$pre_command \n"; }
     if ($opt_v) { print "$regexp \n"; }
     $regexp =~ s/^\s*//;
     $regexp =~ s/\s*$//;
     $lineout = `cd $workdir; $pre_command`;


     if ( $lineout =~ /$regexp/ ) {
	 print "$name: \t [ \033[32m  OK  \033[0m ] \n"; 
         $test_succeded = 1;
     } else {
	 print "$name: \t [ \033[31m FAIL \033[0m ] \n";
         $test_succeded = 0;
     }
   }

 } else {
   if ( $_ =~ /^RUN/) { print " skipping test\n"; }
 }
}


if (!$opt_i) { print "\n\n\n"; }
if ($opt_l && $valid)  { system ("cp $workdir/out out.log"); }
if (!$opt_p && !$opt_m && $test_succeded) { system ("rm -rf $workdir"); }
