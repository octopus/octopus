<html>
<head>
<title>Variables used in octopus</title>
<script>
function op() { //This function is used with folders that do not open pages themselves. See online docs.
}
</script>
</head>

<!--
(Please keep all copyright notices.)
This frameset document includes the Treeview script.
Script found in: http://www.treeview.net
Author: Marcelino Alves Martins

You may make other changes, see online instructions, 
but do not change the names of the frames (treeframe and basefrm)
-->


<FRAMESET cols="200,*" onResize="if (navigator.family == 'nn4') window.location.reload()"> 
  <FRAME src="varsLeftFrame.php?page=<?php echo $_GET["page"]; ?>" name="treeframe" > 
  <FRAME SRC="<?php
if($_GET["section"])
  echo "vars/".$_GET["section"].".html#".$_GET["name"];
else
  echo "varsRightFrame.html"
?>" name="basefrm"> 
</FRAMESET> 


</HTML>
