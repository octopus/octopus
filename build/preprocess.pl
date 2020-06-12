#!/usr/bin/env perl

my $fin         = $ARGV[0];
my $ndebug      = $ARGV[1] eq "no";
my $nline_num   = $ARGV[2] eq "no";

open(IN,  "<$fin");

while($_ = <IN>)
{
  # First, eliminate the push_sub and pop_sub if the code is compiled in non-debug mode.
  # We replace it by an empty line as not to mess up the line numbers.
  if($ndebug){
      if (/push_sub/) {
          print "\n";
          next;
      }
      if (/pop_sub/) {
          print "\n";
          next;
      }
  }

  # Substitute "\newline" by a real new line.
  s/\\newline/\n/g;

  # Substitute "\cardinal" by "#".
  s/\\cardinal/#/g;

  if($nline_num){
      next if /^#/;
  }
  print;
}

close(IN);
