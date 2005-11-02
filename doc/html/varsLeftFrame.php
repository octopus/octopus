<!--
     (Please keep all copyright notices.)
     This frameset document includes the Treeview script.
     Script found at: http://www.treeview.net
     Author: Marcelino Alves Martins

     Instructions:
	 - Through the style tag you can change the colors and types
	   of fonts to the particular needs of your site. 
	 - A predefined block has been made for stylish people with
	   black backgrounds.
-->


<html>

<head>

<style>
   BODY {background-color: white}
   TD {font-size: 10pt; 
       font-family: verdana,helvetica; 
	   text-decoration: none;
	   white-space:nowrap;}
   A  {text-decoration: none;
       color: black}
</style>

<script src="js/ua.js"></script>
<script src="js/ftiens4.js"></script>

<?php
if($_GET["page"] == "alpha")
  echo "<script src='vars/alpha.js'></script>";
else
  echo "<script src='vars/sections.js'></script>";
?>

</head>

<body topmargin=16 marginheight=16>
<div style="text-align:center">
     <?php
if($_GET["page"] == "alpha")
  echo "<a href='vars.php?page=sections' target=_top style='color:blue'>sections</a> | alpha";
else
  echo "sections | <a href='vars.php?page=alpha' target=_top style='color:blue'>alpha</a>";
?>
</div>

<!-- By making any changes to this code you are violating your user agreement.
     Corporate users or any others that want to remove the link should check 
	 the online FAQ for instructions on how to obtain a version without the link -->
<!-- Removing this link will make the script stop from working -->
<div style="position:absolute; top:0; left:0; "><table border=0><tr><td><font size=-2><a style="font-size:7pt;text-decoration:none;color:silver" href="http://www.treemenu.net/" target=_blank>JavaScript Tree Menu</a></font></td></tr></table></div>

<!-- Build the browser's objects and display default view of the 
     tree. -->
<script>initializeDocument()</script>
<noscript>
A tree for site navigation will open here if you enable JavaScript in your browser.
</noscript>

</html>
