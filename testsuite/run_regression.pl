#!/usr/bin/perl
#
# $Id$

use Getopt::Std;


sub usage {

  print <<EndOfUsage;

 Copyright (C) 2005 by Heiko Appel

Usage: run_regression.pl [options]

    -n        dry-run
    -v        verbose
    -h        this usage
    -e        name of octopus executable
    -c        create template
    -f        filename of testsuite
    -d        working directory for the tests
    -i        print inputfile
    -p        preserve working directories

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
Test    :
Author  : $author
Date    : $date
Arch    : $arch
Release :
Enabled : Yes

INP
CalculationMode = gs
fromScratch = yes
%Coordinates
  "Na" | 0 | 0 | -1.54
  "Na" | 0 | 0 |  1.54
%
EOF

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

getopts("nvhe:c:f:d:ip");

# Default values
$workdir = `mktemp -d`;
chomp($workdir);

# Handle options
$opt_h && usage;
$opt_c && create_template;
if($opt_d)  { $workdir = $opt_d; }
if(!$opt_i) { print "Using workdir    : $workdir \n" };

if($opt_e)  { $octopus_exe = $opt_e; }
else{
  $octopus_exe =        "octopus";
  $octopus_exe = "../src/octopus"       if( -x "../src/octopus");
  $octopus_exe = "../src/octopus_debug" if( -x "../src/octopus_debug");
}
if($octopus_exe !~ /^\//){
  $octopus_exe = $ENV{PWD}."/$octopus_exe";
}

die "could not find executable: $octopus_exe \n" if (`which $octopus_exe` eq "");
if(!$opt_i) { print "Using executable : $octopus_exe \n"; }

system ("rm -rf $workdir");
mkdir $workdir;

open(TESTSUITE, "<".$opt_f ) or die "cannot open testsuite file \n";


while ($_ = <TESTSUITE>) {

 # skip comments
 next if /^#/;

 if ( $_ =~ /Test\s*:\s*(.*)\s*$/) {
  %test = ();
  $test{"name"} = $1;
  if(!$opt_i) { print "\033[34m $test{\"name\"} \033[0m \n"; }
 }

 if ( $_ =~ /Author\s*:\s*(.*)\s*$/) {
  %test = ();
  $test{"author"} = $1;
 }

 if ( $_ =~ /Date\s*:\s*(.*)\s*$/) {
  %test = ();
  $test{"date"} = $1;
 }

 if ( $_ =~ /Release\s*:\s*(.*)\s*$/) {
  %test = ();
  $test{"release"} = $1;
 }

 if ( $_ =~ /Enabled\s*:\s*(.*)\s*$/) {
  %test = ();
  $enabled = $1;
  $enabled =~ s/^\s*//;
  $enabled =~ s/\s*$//;
  $test{"enabled"} = $enabled;
 }

 # Running this regression test if it is enabled
 if ( $enabled eq "Yes" ) {

   # generating input file for test run
   if ( $_ =~ /^INP/) {
     $test{"inp"} = 1;
     open(INP,">$workdir/inp");
     while ( ($_ = <TESTSUITE>) ) {
       last if ( $_ =~ /^EOF/ );
       print INP $_;
       if ($opt_i) { print STDOUT $_; }

     }
     close(INP);
   }

   if ( $_ =~ /^RUN/) {
     die "inp not defined" if !$test{"inp"};
     if (!$opt_n) {
       print "\nStarting test run ...\n";
       system("cd $workdir; $octopus_exe < inp &> out");
       print "Finished test run.\n"; }
     else {
       if(!$opt_i) { print "cd $workdir; $octopus_exe < inp &> out \n"; }
     }
     $test{"run"} = 1;
   }

   if ( $_ =~ /^Match/ && !$opt_n) {
     ($match, $name, $pre_command, $regexp) = split(/\;/,$_);

     die "have to run before matching" if !$test{"run"};

     if ($opt_v) { print "$pre_command \n"; }
     if ($opt_v) { print "$regexp \n"; }
     $regexp =~ s/^\s*//;
     $regexp =~ s/\s*$//;
     $lineout = `cd $workdir; $pre_command`;


     if ( $lineout =~ /$regexp/ ) {
       print "$name: \t [ \033[32m  OK  \033[0m ] \n";
     } else {
       print "$name: \t [ \033[31m FAIL \033[0m ] \n";
     }
   }

 } else {
   if ( $_ =~ /^RUN/) { print " skipping test\n"; }
 }
}


if (!$opt_i) { print "\n\n\n"; }
if (!$opt_p) { system ("rm -rf $workdir"); }
