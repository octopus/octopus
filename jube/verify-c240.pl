#!/usr/bin/env perl

my $xmloutfile = $ARGV[0];
my $subdir     = $ARGV[1];
my $vcheck     = 1; # this is the error condition
my $vcomment   = "";

######################################################
# TD stuff

$fname = "$subdir/gs_static/info";

if(not -f $fname){
  $vcomment = "No output file found";
  goto output_xml;
}

$e480ref = -0.202221;
$e480 = `n=\$(cat -n $fname | grep '480   --' | head -n 1 | awk '{print \$1;}'); cat $fname | tail -n +\$((n+0)) |  head -n 1 | cut -b 13- | perl -ne '/\\s*([0-9\-+.eEdD]*)/; print \$1'`;

$vcomment .= "Eigenvalue wrong; " if(abs($e480-$e480ref) > 0.0001);

$etotref = -1372.98162539;
$etot = `n=\$(cat -n $fname | grep 'Total       =' | head -n 1 | awk '{print \$1;}'); cat $fname | tail -n +\$((n+0)) |  head -n 1 | cut -b 20- | perl -ne '/\\s*([0-9\-+.eEdD]*)/; print \$1'`;

$vcomment .= "Total energy wrong; " if(abs($etot-$etotref) > 0.0001);

######################################################
# TD stuff

$fname = "$subdir/td_td.general/energy";

if(not -f $fname){
  $vcomment = "No output file found";
  $vcheck = 1;
  goto output_xml;
}

$tderef = -1.372970304837e+03;
$tde=`cat $fname | tail -n -1 | head -n 1 | cut -b 30- | perl -ne '/\s*([0-9\-+.eEdD]*)/; print \$1'`;

$vcomment .= "Propagation wrong; " if(abs($tde-$tderef) > 0.0001);

######################################################
# finalize and output

if($vcomment eq ""){
  $vcheck   = 0;
  $vcomment = "Result verified";
}

output_xml:
open(XMLOUT,"> $xmloutfile") || die "cannot open file $xmloutfile";
print XMLOUT "<verify>\n";
print XMLOUT " <parm name=\"vcheck\"   value=\"$vcheck\"   type=\"bool\" unit=\"\" />\n";
print XMLOUT " <parm name=\"vcomment\" value=\"$vcomment\" type=\"string\" unit=\"\"/>\n";
print XMLOUT " <parm name=\"ve480\"    value=\"$e480\"     type=\"float\" unit=\"\"/>\n";
print XMLOUT " <parm name=\"ve480ref\" value=\"$e480ref\"  type=\"float\" unit=\"\"/>\n";
print XMLOUT " <parm name=\"vetot\"    value=\"$etot\"     type=\"float\" unit=\"\"/>\n";
print XMLOUT " <parm name=\"vetotref\" value=\"$etotref\"  type=\"float\" unit=\"\"/>\n";
print XMLOUT " <parm name=\"vtde\"     value=\"$tde\"      type=\"float\" unit=\"\"/>\n";
print XMLOUT " <parm name=\"vtderef\"  value=\"$tderef\"   type=\"float\" unit=\"\"/>\n";
print XMLOUT "</verify>\n";
print XMLOUT "\n";
close(XMLOUT);

exit($vcheck);
