---
title: "Installation"
series: "Manual"
weight: 4
---


Maybe somebody else installed {{< octopus >}} for you. In that case, the files
should be under some directory that we can call {{< inst_file >}}, the executables
in {{< inst_file "bin/" >}} (e.g. if {{< inst_file >}}={{< file "/usr/local/" >}}, the main {{< octopus >}}
executable is then in {{< file "/usr/local/bin/octopus" >}}); the pseudopotential files that {{< octopus >}} will
need in {{< inst_file "share/octopus/PP/" >}}, etc.

However, you may be unlucky and that is not the case. In the following we will try
to help you with the still rather unfriendly task of compiling and installing the {{< octopus >}}.

###  Instructions for Specific Architectures  

See {{< Manual "installation/Specific architectures" "step-by-step_instructions" >}} for some specific supercomputers and generic configurations (including Ubuntu and Mac OSX). If the system you are using is in the list, this will be the easiest way to install. Add entries for other supercomputers or generic configurations on which you have successfully built the code.

##  Dowloading  

Download the latest {{< octopus >}} version here: [[Octopus {{< Octopus major version >}}]].

###  Binaries  

If you want to install {{< octopus >}} in a {{< name "Debian" >}}-based {{< name "Linux" >}} box ({{< name "Debian" >}} or {{< name "Ubuntu" >}}, you might not need to compile it;  we release binary packages for some platforms. Keep in mind that these packages are intended to run on different systems and are therefore only moderately optimized. If this is an issue you have to compile the package yourself with the appropriate compiler flags and libraries (see below).

Download the appropriate {{< file ".deb" >}} file from the [[Octopus {{< octopus major version >}} | downloads page for the current version]]. Install it (using root access) with the command below, using the appropriate filename you downloaded:

<!--- {{< command_line "dpkg -i octopus_package.deb" "<nowiki>-</nowiki>" >}} -->
```bash
dpkg -i octopus_package.deb
```

###  Source code  

If you have a different system from those mentioned above or you want to compile {{< octopus >}} you need to get the source code file ({{< file ".tar.gz" >}} file) and follow the compilation instructions below.

##  Building  
###  Quick instructions  

For the impatient, here is the quick-start: 
<!-- 
{{< command_line "tar xzf octopus-{{< octopus_version >" >}}.tar.gz}} 
{{< command_line "cd octopus-{{< octopus_version >" >}}}}
{{< command_line "./configure" >}}
{{< command_line "make" >}}
{{< command_line "make install" >}}
-->

```bash
tar xzf octopus-{{< octopus_version >}}.tar.gz
cd octopus-{{< octopus_version >}}
./configure
make
make install
```

This will probably {{< emph "not" >}} work, so before giving up, just read
the following paragraphs.

###  Slow instructions  

There is an {{< Manual "installation/Building from scratch" "appendix" >}} with detailed instructions on how to compile {{< octopus >}} and the required libraries from scratch -- you only need to do this if you are unable to install from a package manager or use per-built libraries.

###  Long instructions  

The code is written in standard {{< name "Fortran 2003" >}}, with some routines written in {{< name "C" >}} (and in {{< name "bison" >}}, if we count the input parser). To build it you will need both a {{< name "C" >}} compiler ({{< name "gcc" >}} works just fine and it is available for almost every piece of silicon), and a
{{< name "Fortran 2003" >}} compiler. You can check in the {{< Manual "Appendix/Compilers" "Compilers Appendix" >}} which compilers {{< octopus >}} has been tested with. This appendix also contains hints on potential problems with certain platform/compiler combinations and how to fix them.

####  Requirements  
<br>

Besides the compiler, you will also need:

* {{< name "make" >}}:  most computers have it installed, otherwise just grab and install the {{< name "GNU make" >}}.

* {{< name "cpp" >}}: The {{< name "C" >}} preprocessor is heavily used in {{< octopus >}} to preprocess {{< name "Fortran" >}} code. It is used for both C (from the CPP variable) and Fortran (FCCPP). {{< name "GNU cpp" >}} is the most convenient but others may work too. For more info, see [[Preprocessors]].

* {{< name "Libxc" >}}: The library of exchange and correlation functionals. It used to be a part of {{< octopus >}}, but since version 4.0.0 it is a standalone library and needs to be installed independently.  For more information, see the [http://www.tddft.org/programs/Libxc Libxc page]. {{< octopus >}} 4.0.0 and 4.0.1 require version 1.1.0 (not 1.2.0 or 1.0.0). {{< octopus >}} 4.1.2 requires version 2.0.x or 2.1.x, and won't compile with 2.2.x. (Due to bugfixes from libxc version 2.0 to 2.1, there will be small discrepancies in the testsuite for {{< code "functionals/03-xc.gga_x_pbea.inp" >}} and {{< code "periodic_systems/07-tb09.test" >}}). {{< octopus >}} 5.0.0 supports libxc versions 2.0.x, 2.1.x and 2.2.x. Please note: The Libxc testsuite prior to 2.1 will report some errors in most cases. This is not something to worry about.

* {{< name "FFTW" >}}: We have relied on this great library to perform Fast Fourier Transforms ({{< name "FFTs" >}}). You may grab it from the [http://www.fftw.org/ {{< name "FFTW" >}} site]. You require {{< name "FFTW" >}} version 3.

* {{< name "LAPACK/BLAS" >}}: Our policy is to rely on these two libraries as much as possible on these libraries for linear-algebra operations. If you are running {{< name "Linux" >}}, there is a fair chance they are already installed in your system. The same goes to the more heavy-weight machines ({{< name "alphas" >}}, {{< name "IBMs" >}}, {{< name "SGIs" >}}, etc.). Otherwise, just grab the source from [http://www.netlib.org {{< name "netlib" >}} site].

* {{< name "GSL" >}}: Finally someone had the nice idea of making a public scientific library! {{< name "GSL" >}} still needs to grow, but it is already quite useful and impressive. {{< octopus >}} uses splines, complex numbers, special functions, etc. from {{< name "GSL" >}}, so it is a must! If you don't have it already installed in your system, you can obtain {{< name "GSL" >}} from the [http://www.gnu.org/software/gsl/ {{< name "GSL" >}} site]. You will need version 1.9 or higher. Version 4.0 of {{< octopus >}} (and earlier) can only use GSL 1.14 (and earlier). A few tests will fail if you use GSL 1.15 or later. Version 5.0.0 of {{< octopus >}} (and earlier) can only use GSL 1.16 or earlier, due to a bug in our configure script.

* {{< name "Perl" >}}: During the build process {{< octopus >}} runs several scripts in this language. It's normally available in every modern {{< name "Unix" >}} system.
<br>

####  Optional libraries  
<br>

There are also some optional packages; without them some parts of {{< octopus >}} won't work:
<br>
* {{< name "MPI" >}}: If you want to run {{< octopus >}} in multi-tentacle (parallel) mode, you will need an implementation of {{< name "MPI" >}}. [http://www-unix.mcs.anl.gov/mpi/mpich/ MPICH] or [http://www.open-mpi.org/ Open MPI] work just fine in our {{< name "Linux" >}} boxes.

* {{< name "PFFT" >}}: We rely on this great library for highly scalable parallel Poisson solver, based on Fast Fourier Transforms ({{< name "FFTs" >}}). You may grab it from the [https://www-user.tu-chemnitz.de/~potts/workgroup/pippig/software.php.en {{< name "M. Pippig's" >}} site]. You also require {{< name "FFTW" >}} version 3.3 compiled with MPI and with a small patch by M. Pippig (also available there).

* {{< name "NetCDF" >}}: The [http://www.unidata.ucar.edu/software/netcdf/ Network Common Dataform] library is needed for writing the binary files in a machine-independent, well-defined format, which can also be read by [[Manual:Visualization|visualization programs]] such as [http://www.opendx.org/ OpenDX]

* {{< name "GDLib" >}}: A library to read graphic files. See [[Tutorial:Particle in an octopus]]. (The simulation box in 2D can be specified via {{< Variable2 "BoxShapeImage" >}}.) Available from http://www.libgd.org/

* {{< name "SPARSKIT" >}}: [http://www-users.cs.umn.edu/~saad/software/SPARSKIT/ Library for sparse matrix calculations]. Used for one propagator technique.

* {{< name "ETSF I/O" >}}: An input/output library implementing the ETSF standardized formats, requiring NetCDF, available at [http://www.etsf.eu/resources/software/libraries_and_tools]. Versions 1.0.2, 1.0.3, and 1.0.4 are compatible with {{< octopus >}} (though 1.0.2 will produce a small discrepancy in a filesize in the testsuite). It must have been compiled with the same compiler you are using with {{< octopus >}}. To use ETSF_IO, include this in the <tt>configure</tt> line, where <tt>$DIR</tt> is the path where the library was installed:

```text
 --with-etsf-io-prefix="$DIR"
```

* {{< name "LibISF" >}}: (version 5.0.0 and later) To perform highly scalable parallel Poisson solver, based on BigDFT 1.7.6, with a cheap memory footprint. You may grab it from the [http://bigdft.org/ {{< name "BigDFT" >}} site]. You require {{< name "BigDFT" >}} version 1.7.6 compiled with MPI, following these instructions: [http://bigdft.org/Wiki/index.php?title=Installation-Building_the_Poisson_Solver_library_only installation instructions ]. Probably, you have to manually copy the files "libwrappers.a" and "libflib.a" to the installation "/lib" directory. To configure {{< octopus >}}, you have to add this configure line:

```text
  --with-isf-prefix="$DIR"
```

####  Unpacking the sources  

Uncompress and untar it ({{< command "gzip -cd octopus-{{< octopus_version >" >}}.tar.gz | tar -xvf -}}). In the following, {{< file "OCTOPUS-SRC/" >}} denotes the source directory of {{< octopus >}}, created by the {{< command "tar" >}} command.


The {{< file "OCTOPUS-SRC/" >}} contains the following subdirectories of interest to users:

;{{< file "doc/" >}}: The documentation of {{< octopus >}}, mainly in HTML format.

;{{< file "liboct_parser/" >}}: The C library that handles the input parsing.

;{{< file "share/PP/" >}}: Pseudopotentials. In practice now it contains the Troullier-Martins (PSF and UPF formats) and Hartwigsen-Goedecker-Hutter pseudopotential files.

;{{< file "share/util/" >}}: Currently, the {{< emph "utilities" >}} include a couple of IBM OpenDX networks ({{< file "mf.net" >}}), to visualize wavefunctions, densities, etc.

;{{< file "testsuite/" >}}: Used to check your build. You may also use the files in here as samples of how to do various types of calculations.

;{{< file "src/" >}}: Fortran90 and C source files. Note that the Fortran90 files have to be preprocessed before being fed to the Fortran compiler, so do not be scared by all the - directives.

####  Development version  

You can get the development version of {{< octopus >}} by downloading it from the {{< octopus >}} project on [https://gitlab.com/octopus-code/octopus gitlab.com].

You can also get the current version with the following command (you need the {{< name "git" >}} package):

{{< command_line "git clone git@gitlab.com:octopus-code/octopus.git" >}}

Before running the configure script, you will need to run the GNU autotools. This may be done by executing:

{{< command_line "autoreconf -i" >}}

Note that you need to have working recent versions of the {{< name "automake" >}} and {{< name "autoconf" >}}. In particular, the configure script may fail in the part <tt>checking for Fortran libraries of mpif90</tt> for <tt>autoconf</tt> version 2.59 or earlier. The solution is to update <tt>autoconf</tt> to 2.60 or later, or manually set <tt>FCLIBS</tt> in the <tt>configure</tt> command line to remove a spurious apostrophe.

If autoreconf is failing with "aclocal: warning: couldn't open directory 'm4': No such file or directory", create an empty folder named m4 inside external_libs/spglib-1.9.9/.

Please be aware that the development version may contain untested changes that can affect the execution and the results of {{< octopus >}}, especially if you are using new and previously unreleased features. So if you want to use the development version for production runs, you should at least contact {{< octopus >}} developers.

####  Configuring  

Before configuring you can (should) set up a couple of options. Although
the {{< name "configure" >}} script tries to guess your system settings for you, we recommend 
that you set explicitly the default Fortran compiler and the compiler options.
Note that {{< name "configure" >}} is a standard tool for Unix-style programs and you can find a lot of generic documentation on how it works elsewhere.

For example, in {{< name "bash" >}} you would typically do:

{{< command_line "export FC<nowiki>=</nowiki>ifort" >}}
{{< command_line "export FCFLAGS<nowiki>=</nowiki>&quot;-O2 -xHost&quot;" >}}

if you are using the Intel Fortran compiler on a linux machine.

Also, if you have some of the required libraries in some unusual directories,
these directories may be placed in the variable {{< code "LDFLAGS" >}} (e.g.,
{{< command "export LDFLAGS<nowiki>=</nowiki>$LDFLAGS:/opt/lib/" >}}).

The configuration script will try to find out which compiler you are using.
Unfortunately, and due to the nature of the primitive language that {{< octopus >}}
is programmed in, the automatic test fails very often. Often it is better to set
the variable {{< code "FCFLAGS" >}} by hand, check the [[Manual:Appendix:Compilers|Compilers Appendix]] page for which flags have been reported to work with different Fortran compilers.

You can now run the configure script 
{{< command_line "./configure" >}}

You can use a fair amount of options to spice {{< octopus >}} to your own taste. To obtain a full list just type {{< command "./configure --help" >}}. Some 
commonly used options include:

;{{< flag "<nowiki>--prefix=</nowiki>" >}}{{< inst_file >}}: Change the base installation dir of {{< octopus >}} to {{inst_file|}}. {{< inst_file >}} defaults to the home directory of the user who runs the {{< name "configure" >}} script.

;{{< flag "--with-fft-lib<nowiki>=</nowiki>" >}}{{< file "&lt;lib&gt;" >}}: Instruct the {{< name "configure" >}} script to look for the {{< name "FFTW" >}} library exactly in the way that it is specified in the {{< file "&lt;lib&gt;" >}} argument. You can also use the {{< code "FFT_LIBS" >}} environment variable.

;{{< flag "--with-pfft-prefix<nowiki>=</nowiki>" >}}{{< file "DIR/" >}}: Installation directory of the {{< name "PFFT" >}} library.

;{{< flag "--with-pfft-lib<nowiki>=</nowiki>" >}}{{< file "&lt;lib&gt;" >}}: Instruct the {{< name "configure" >}} script to look for the {{< name "PFFT" >}} library exactly in the way that it is specified in the {{< file "&lt;lib&gt;" >}} argument. You can also use the {{< code "PFFT_LIBS" >}} environment variable.

;{{< flag "--with-blas<nowiki>=</nowiki>" >}}{{< file "&lt;lib&gt;" >}}: Instruct the {{< name "configure" >}} script to look for the {{< name "BLAS" >}} library in the way that it is specified in the {{< file "&lt;lib&gt;" >}} argument.

;{{< flag "--with-lapack<nowiki>=</nowiki>" >}}{{< file "&lt;lib&gt;" >}}: Instruct the {{< name "configure" >}} script to look for the {{< name "LAPACK" >}} library in the way that it is specified in the {{< file "&lt;lib&gt;" >}} argument.

;{{< flag "--with-gsl-prefix<nowiki>=</nowiki>" >}}{{< file "DIR/" >}}: Installation directory of the {{< name "GSL" >}} library. The libraries are expected to be in {{< file "DIR/lib/" >}} and the include files in {{< file "DIR/include/" >}}. The value of {{< file "DIR/" >}} is usually found by issuing the command {{< command "gsl-config --prefix" >}}.

;{{< flag "--with-libxc-prefix<nowiki>=</nowiki>" >}}{{< file "DIR/" >}}: Installation directory of the {{< name "Libxc" >}} library.

If you have problems when the {{< name "configure" >}} script runs, you can find more details of what happened in the file {{< name "config.log" >}} in the same directory.

####  Compiling and installing  

Run {{< command "make" >}} and then {{< command "make install" >}}. The compilation may take some time, so you might want to speed it up by running {{< command "make" >}} in parallel ({{< command "make -j" >}}). If everything went fine, you should now be able to taste {{< octopus >}}. 

Depending on the value given to the {{< name "--prefix" >}}={{< inst_file >}} given, the executables will reside in {{< inst_file "bin/" >}}, and the auxiliary files will be copied to {{< inst_file "share/octopus" >}}.

####  Testing your build  

After you have successfully built {{< octopus >}}, to check that your build works as expected there is a battery of tests that you can run.  They will check that {{< octopus >}} executes correctly and gives the expected results (at least for these test cases). If the parallel version was built, the tests will use up to 6 MPI processes, though it should be fine to run on only 4 cores. (MPI implementations generally permit using more tasks than actual cores, and running tests this way makes it likely for developers to find race conditions.)

To run the tests, in the sources directory of {{< octopus >}} use the command

{{< command_line "make check" >}}

or if you are impatient,

{{< command_line "make check-short" >}}

which will start running the tests, informing you whether the tests are passed or not. For examples of job scripts to run on a machine with a scheduler, please see {{< Manual "Specific_architectures" "Specific_architectures" >}}.

If all tests fail, maybe there is a problem with your executable (like a missing shared library). 

If only some of the tests fail, it might be a problem when calling some external libraries (typically blas/lapack).  Normally it is necessary to compile all Fortran libraries with the same compiler. If you have trouble, try to look for help in the [http://www.tddft.org/mailman/listinfo/octopus-users  {{< octopus >}} mailing list].

####  Fast recompilation  

NOTE: This feature is currently only available in the development version and the plan is to include it in the {{< Octopus >}} 9 release.

If you have already compiled the code and if you are changing only one file, you can run

{{< command_line "make NODEP<nowiki>=</nowiki>1" >}}

to ignore the dependencies due to Fortran module files. This will only compile the files that have changed and link the executables; therefore, it is much faster. If you change, e.g., interfaces of modules or functions, you need to to run {{< command "make" >}} without {{< command "NODEP<nowiki>=</nowiki>1" >}} to ensure a correct handling of the dependencies.


{{< manual_foot prev="Manual:About Octopus" next="Manual:Input file" >}}
---------------------------------------------
