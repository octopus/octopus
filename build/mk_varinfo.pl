#!/usr/bin/perl

if(!-d 'src' && !-d 'share'){
  print stderr "Please run this script from the octopus main directory\n";
  exit;
}

opendir(DIR, 'src');
@F90 = grep { /\.F90$/ && -f "src/$_" } readdir(DIR);
closedir DIR;

open(OUT, ">share/varinfo");

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

	if(/^ /){ # get whole description
	  $desc .= " ".$_;
	}else{
	  if($desc ne "") {
	    print_desc($desc);
	    $desc = "";
	  }

	  printf OUT "%s\n", $_;
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
	printf OUT "END\n\n";
      }
    }
  }
  close(IN);
}

close(OUT);

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

  $desc =~ s/^\s*//;
  $desc =~ s/\s\s+/ /g;

  while(){
    if(length($desc) <= $ml){
      print OUT " ", $desc, "\n";
      last;
    }

    for($i=$ml; $i>0 && substr($desc, $i, 1) ne ' '; $i--){};
    $line = substr($desc, 0, $i);
    $desc = substr($desc, $i+1);

    $spaces = $ml - length($line);
    for($i=0; $i<$spaces; $i++){
      $line =~ s/([^ ]) ([^ ])/$1  $2/;
    }
    print OUT " ", $line, "\n";
  }
}
