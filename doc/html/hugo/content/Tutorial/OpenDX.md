---
title: "OpenDX"
tags: ["Obsolete"]
series: "Tutorial"
---


''Note: the OpenDX program is no longer supported and the following tutorial should be considered obsolete. We recommend you to use one of the alternatives mentioned in this [tutorial](../Benzene_molecule).''

##  Input  

At this point, the input should be quite familiar to you:

```text
 {{< Variable2 "CalculationMode" >}} = gs
 {{< Variable2 "UnitsOutput" >}} = eV_Angstrom
 
 {{< Variable2 "Radius" >}} = 5*angstrom
 {{< Variable2 "Spacing" >}} = 0.15*angstrom
 
 {{< Variable2 "Output" >}} = wfs + density + elf + potential
 {{< Variable2 "OutputFormat" >}} = cube + xcrysden + dx + axis_x + plane_z
 
 {{< Variable2 "XYZCoordinates" >}} = "benzene.xyz"
```

Coordinates are in this case given in the file {{< file "benzene.xyz" >}} file. This file should look like:

```text
 12
    Geometry of benzene (in Angstrom)
 C  0.000  1.396  0.000
 C  1.209  0.698  0.000
 C  1.209 -0.698  0.000
 C  0.000 -1.396  0.000
 C -1.209 -0.698  0.000
 C -1.209  0.698  0.000
 H  0.000  2.479  0.000
 H  2.147  1.240  0.000
 H  2.147 -1.240  0.000
 H  0.000 -2.479  0.000
 H -2.147 -1.240  0.000
 H -2.147  1.240  0.000
```

Note that we asked octopus to output the wavefunctions (<tt>wfs</tt>), the density, the electron localization function (<tt>elf</tt>) and the Kohn-Sham potential. We ask Octopus to generate this output in several formats, that are requires for the different visualization tools that we mention below. If you want to save some disk space you can just keep the option that corresponds to the program you will use. Take a look at the documentation on variables {{< Variable2 "Output" >}} and {{< Variable2 "OutputFormat" >}} for the full list of possible quantities to output and formats to visualize them in.

###  OpenDX  

[http://www.opendx.org OpenDX] is a very powerful, although quite complex, program for the 3D visualization of data. It is completely general-purpose, and uses a visual language to process the data to be plotted. Don't worry, as you will probably not need to learn this language. We have provided a program ({{< file "mf.net" >}} and {{< file "mf.cfg" >}}) that you can get from the {{< file "SHARE/util" >}} directory in the {{< octopus >}} installation (probably it is located in {{< file "/usr/share/octopus/util" >}}). 

Before starting, please make sure you have it installed together with the Chemistry extensions. You may find some instructions [[Releases-OpenDX|here]]. Note that openDX is far from trivial to install and to use. However, the outcomes are really spectacular, so if you are willing to suffer a little at the beginning, we are sure that your efforts will be compensated.

If you run {{< octopus >}} with this input file, afterward you will find in the {{< file "static" >}} directory several files with the extension {{< file ".dx" >}}. These files contain the wave-functions, the density, etc. in a format that can be understood by [http://www.opendx.org/ Open DX] (you will need the chemistry extensions package installed on top of DX, you can find it in our [[releases]] page). 

So, let us start. First copy the files {{< file "mf.net" >}} and {{< file "mf.cfg" >}} to your open directory, and open dx with the command line

```text
  > dx
```

If you followed the instructions given in [[Releases-OpenDX|here]] you should see a window similar to the one shown on the right.
<gallery>
```text
  Image:Tutorial_dx_1.png|Screenshot of the Open DX main menu
  Image:Tutorial_dx_2.png|Main dialog of mf.net
  Image:Tutorial_dx_3.png|Isosurfaces dialog
  Image:Tutorial_dx_4.png|Geometry dialog
```
</gallery>

Now select '''<tt>Run Visual Programs...</tt>''', select the {{< file "mf.net" >}} file, and click on '''<tt>OK</tt>'''. This will open the main dialog of our application. You will no longer need the main DX window, so it is perhaps better to minimize it to get it out of the way. Before continuing, go to the menu '''<tt>Execute > Execute on Change</tt>''', which instructs to open DX to recalculate out plot whenever we make a change in any parameters.

{{< figure src="/images/Tutorial_dx_5.png" width="500px" caption="350px" >}}

Now, let us choose a file to plot. As you can see, the field <tt>File</tt> in the mains dialog is <tt>NULL</tt>. If you click on the three dots next to the field a file selector will open. Just select the file {{< file "static/density-1.dx" >}} and press OK. Now we will instruct open DX to draw a iso-surface. Click on <tt>Isosurfaces</tt> and a new dialog box will open (to close that dialog box, please click again on <tt>Isosurfaces</tt> in the main dialog!) In the end you may have lots of windows lying around, so try to be organized. Then click on <tt>Show Isosurface 1</tt> - an image should appear with the density of benzene. You can zoom, rotate, etc the picture if you follow the menu '''<tt>Options > View control...</tt>''' in the menu of the image panel. You can play also around in the isosurfaces dialog with the <tt>Isosurface value</tt> in order to make the picture more appealing, change the color, the opacity, etc. You can even plot two isosurfaces with different values.

We hope you are not too disappointed with the picture you just plotted. Densities are usually quite boring (even if the Hohenberg-Kohn theorem proves that they know everything)! You can try to plot other quantities like the electronic localization function, the wave-functions, etc, just by choosing the appropriate file in the main dialog. You can try it later: for now we will stick to the "boring density".

Let us now add a structure. Close the isosurfaces dialog '''by clicking again in <tt>Isosurfaces</tt> in the main dialog box''', and open the <tt>Geometry</tt> dialog. In the <tt>File</tt> select the file {{< file "benzene.xyz" >}}, and click on <tt>Show Geometry</tt>. The geometry should appear beneath the isosurface. You can make bonds appear and disappear by using the slide <tt>Bonds length</tt>,  change its colors, etc. in this dialog.

To finish, let us put some "color" in our picture. Close the <tt>Geometry</tt> dialog (you know how), and open again the isosurfaces dialog. Choose <tt>Map isosurface</tt> and in the <tt>File</tt> field choose the file {{< file "static/vh.dx" >}}. There you go, you should have a picture of an isosurface of benzene colored according to the value of the electrostatic potential! Cool, no? If everything went right, you should end with a picture like the one on the right.

As an exercise you can try to plot the wave-functions, etc, make slabs, etc. Play around. As usual some things do not work perfectly, but most of it does!


---------------------------------------------
