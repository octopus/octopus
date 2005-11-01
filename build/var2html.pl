#!/usr/bin/perl

# configuration
$out_dir = "html";

if(! -f 'share/varinfo'){
  print stderr "Please run this script from the octopus main directory\n";
  exit;
}

# let us parse the varinfo file
open(IN, '<share/varinfo');
$thisvar = "";
$thisfield = "";
while($l = <IN>){
  if($thisvar && $thisfield && $l !~ /^\w/){
    if($l =~ /^\s*$/){
      $vars{$thisvar}{$thisfield} .= "\n<br>\n";
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
    $vars{$thisvar}{"sort"}     = $thisvar;
  }

  next if(!$thisvar);

  if($key eq "section"){
    $vars{$thisvar}{"section"} = $arg;
    $vars{$thisvar}{"sort"} = $vars{$thisvar}{"section"}."::".$vars{$thisvar}{"variable"};
  }elsif($key eq "type"){
    $vars{$thisvar}{"type"} = $arg;
  }elsif($key eq "default"){
    $vars{$thisvar}{"default"} = $arg;
  }elsif($key eq "description"){
    $thisfield = "description";
  }elsif($key eq "option"){
    $arg =~ /^(\S*)\s*(.*?)\s*$/;
    $thisfield = "option_$2_$1";
  }elsif($key eq "end"){
    $thisvar   = "";
    $thisfield = "";
  }
}
close(IN);

@sorted = sort { $vars{$a}{"sort"} cmp $vars{$b}{"sort"} } keys %vars;

# make html pages
`rm -rf $out_dir`;
mkdir($out_dir);
print_index();
print_vars();

sub print_index(){
  my $i, $j, $k, $key, $sect, $new_sect;
  my @old, @new;

  open(OUT, ">$out_dir/index.html");

  print OUT "<ul>\n";

  $sect = "";
  for($i=0; $i<=$#sorted; $i++){
    $key = $sorted[$i];

    $new_sect = lc($vars{$key}{"section"});
    if($sect ne $new_sect){
      @old = split("::", $sect);
      @new = split("::", $new_sect);

      for($j=0; $old[$j] && $new[$j] && $old[$j] eq $new[$j]; $j++){};

      for($k=$j; $old[$k]; $k++){
	print OUT "</ul>\n";
      }

      for($k=$j; $new[$k]; $k++){
	print OUT "<li>", $new[$k], "</li>\n<ul>\n";
      }

      $sect = $new_sect;
    }

    $sect =~ /^([^:]*)/;
    $link = $1.".html#".$vars{$key}{"variable"};
    $link =~ s/\s/_/g;

    print OUT "<li><a href='$link'>", $vars{$key}{"variable"}, " </a></li>\n";
  }

  #print last <\ul>
  @old = split("::", $sect);
  for($j=0; $old[$j]; $j++){
    print OUT "</ul>\n";
  }
  print OUT "</ul>\n";

  close(OUT);
}

sub print_vars(){
  my $i, $k, $sect, $new_sect;

  $sect = "";
  for($i=0; $i<=$#sorted; $i++){
    $key = $sorted[$i];

    $new_sect = lc($vars{$key}{"section"});
    $new_sect =~ s/:.*$//;
    $new_sect =~ s/\s/_/g;

    if($sect ne $new_sect){
      my $link;

      close(OUT) if($sect);
      $sect = $new_sect;

      open(OUT, ">$out_dir/$sect.html");
    }

    # print the information about the variable
    print OUT 
      "\n
<p><b><a name='$vars{$key}{variable}'>$vars{$key}{variable}</a></b>
<br/><i>Section</i>: $vars{$key}{section}
<br/><i>Type</i>: $vars{$key}{type}
<br/><i>Default</i>: $vars{$key}{default}
<br/>$vars{$key}{description}
";

    my $first = 1;
    foreach $k (sort keys %{ $vars{$key} }) {
      if($k =~ /^option_(.*?)_(.*)/) {
	if($first){
	  print OUT "<br/><i>Options</i>:\n<ul>\n";
	  $first = 0;
	}
	print OUT "<li><b>$2</b>";
	print OUT " ($1)" if($1);
	print OUT ": ", $vars{$key}{$k}, "</li>\n";
      }
    }
    print OUT "</ul>\n" if(!$first);
    print OUT "</p><hr/>\n";
  }
  close(OUT);
}
