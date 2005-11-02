#!/usr/bin/perl

if(!-d 'src' && !-d 'share'){
  print stderr "Please run this script from the octopus main directory\n";
  exit;
}

opendir(DIR, 'src');
@F90 = grep { /\.F90$/ && -f "src/$_" } readdir(DIR);
closedir DIR;

open(OUT_text, ">share/varinfo");
open(OUT_orig, ">share/varinfo_orig");

%opt = ();
foreach $F90file (@F90){
  open(IN, "<src/$F90file");
  while($_=<IN>){
    if(/!%Variable\s+(\S+)/i){
      $var = $1;

      $desc = "";
      do {
	s/^\s*!%//; s/\s*$//;

	if(/^Option\s+(\S+)\s+(\S+)/){
	  put_opt($1, $2);
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
	exit;
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


#####################################################
# tries to put an option in global %opt
sub put_opt{
  my ($a, $b) = @_;

  if($opt{$a} && ($opt{$a} ne $b)){
    print stderr "Option '", $a, "' is multiply defined\n",
      "    '", $opt{$b}, "' ne '", $b, "'\n";
    exit;
  }
  $opt{$a} = $b;
}


#####################################################
# reads in share/variables.local, and then prints %opt to share/variables
sub print_opt{
  # first read in variables.local file
  open(IN, "<share/variables.local");
  while($_=<IN>){
    if(/^\s*(\S+)\s*=\s*(\S+)/){
      put_opt($1, $2);
    }
  }
  close(IN);

  # now print all variables for octopus
  open(OUT, ">share/variables");
  my $key;
  foreach $key (sort(keys %opt)) {
    print OUT $key, " = ", $opt{"$key"}, "\n";
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
