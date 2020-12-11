---
title: "Building from scratch"
series: "Manual"
---


Compiling {{< octopus >}} on some systems might be difficult, this is mainly because some libraries are required and they have to be all compiled with same Fortran compiler. This is a guide on how to compile {{< octopus >}} and the required libraries starting from scratch (a Fortran compiler is required) in a standard Unix system.

'''NOTE: This page should be considered a last resort. On most systems you will be able to install many of these dependencies from a package manager or other pre-built binaries, which will save you a lot of effort.'''

###  Prerequisites  

You will need the following:

* A standard Unix/Linux system, where you want to install {{< octopus >}}.
* A basic knowledge on how to use the system (manage files, edit text files, extract tar archives, etc.).
* A directory where we will install the created binaries, this can be a specially created directory or even your {{< code "$HOME" >}}. From now, we will assume that the directory is {{< basedir >}}, whenever it appears you should replace it with the actual directory you are using. ''Note that the files that you will download and the directories with the sources, including the ones of {{< octopus >}}, don't need to be placed in this directory.''
* A working C compiler, {{< name "gcc" >}} is the ideal choice for {{< name "x86/x86_64" >}} systems. For other platforms you might want to choose a compiler that produces more optimized code. We will assume that the command to run the compiler is {{< cc >}}.
* A working Fortran 95 compiler, if you don't have one you can probably get [http://gcc.gnu.org/fortran/ {{< name "gfortran" >}}] or [http://www.g95.org/ {{< name "g95" >}}]. We will assume that the command to run the Fortran compiler is {{< f90 >}}.
* If you are compiling for a dual 32/64 bits architecture, it is a good idea to tell the compilers explicitly what type of binaries you want, many problems are caused by unadvertedly mixing 32 and 64 bits code. For {{< octopus >}} normally you will want a 64 bits binary ({{< code "-m64" >}} flag in most compilers). Since this flag should be passed always we recommend you to include it in the compiler command used in this guide (sometimes optimization flags are left out), so for example {{< cc >}} could be {{< command "gcc -m64" >}}.

###  Compiler flags  

Since probably you want {{< octopus >}} to run fast, you will probably would to set optimization flags for both the C and Fortran compilers on your system (at least -O3). In general {{< octopus >}} should run fine with aggressive optimization options. We will assume now that you have found a set of flags for {{< cc >}} and {{< f90 >}}, we will call them {{< cflags >}} and {{< fflags >}}. For example {{< cflags >}} could be {{< command "<nowiki>-O3 -march=native</nowiki>" >}}.

''Note: if {{< cc >}}, {{< f90 >}}, {{< cflags >}} or {{< fflags >}} contain spaces you should not enclose them in "". We will put the "" explicitly when it is necessary.''

###  Compilation of libraries  

First we will compile and install the libraries required by {{< octopus >}}.
####  BLAS  

The first library we will compile is {{< name "blas" >}}. We will use the generic reference implementation, this version is not highly optimized but it is free software and simple to install. So after you have a working version of {{< octopus >}} you may want to use a more optimized {{< name "blas" >}} implementation as [http://www.tacc.utexas.edu/resources/software/-blas {{< name "libgoto" >}}], [http://math-atlas.sourceforge.net {{< name "ATLAS" >}}], {{< name "MKL" >}}, {{< name "ACML" >}}, {{< name "ESSL" >}}, etc.

* Download the package from the {{< name "netlib" >}} site: http://www.netlib.org/blas/blas.tgz
* Extract the package and enter into the newly created {{< file "BLAS" >}} directory.
* Edit the {{< file "make.inc" >}} file and modify the definition of the fortran compiler that now should be:

```text
 FORTRAN  = {{< f90 >}}
 OPTS     = {{< fflags >}}
 DRVOPTS  = $(OPTS)
 NOOPT    =
 LOADER   = {{< f90 >}}
 LOADOPTS =
```
* Now type {{< command "make" >}} to compile, this will take a while.
* One compilation is finished, create the directory {{< basedir >}}{{< file "/lib" >}} 
```text
 mkdir {{< basedir >}}/lib
```
and copy the file {{< file "blas_LINUX.a" >}} to {{< basedir >}}{{< file "/lib/libblas.a" >}}
```text
 cp blas_LINUX.a {{< basedir >}}/lib/libblas.a
```

''Note: the generated file librrary will be always called {{< file "blas_LINUX.a" >}} independently of the operating system you are using, this is just a name.''

####  LAPACK  

We will use the open source reference implementation of Lapack.

* Get the Lapak package from netlib: http://www.netlib.org/lapack/lapack.tgz
* Extract the archive and enter the newly created lapack directory.
* Copy {{< file "make.inc.example" >}} to {{< file "make.inc" >}}:
```text
 cp make.inc.example make.inc
```
* Edit {{< file "make.inc" >}} to indicate the compiler the will be used. The relevant part of the file should look like:
```text
 FORTRAN  = {{< f90 >}}
 OPTS     = {{< fflags >}}
 DRVOPTS  = $(OPTS)
 NOOPT    =
 LOADER   = {{< f90 >}}
 LOADOPTS =
```
* Now build the library with
```text
 make lib
```
* Copy the newly created library {{< file "lapack_LINUX.a" >}} to {{< file "<basedir>/lib/liblapack.a" >}} .

####  GSL  

* Get the source of the latest version of GSL from ftp://ftp.gnu.org/gnu/gsl/ , currently it is GSL 1.16: ftp://ftp.gnu.org/gnu/gsl/gsl-1.16.tar.gz .
* Extract the archive and enter the newly created directory.
* Run the configure script with:
```text
 ./configure CC="{{< cc >}}" --prefix={{< basedir >}} --disable-shared --enable-static
```
* Compile and install:
```text
 make
 make install
```

####  FFTW 3  

* Get the sources for the latest version of [http://www.fftw.org/ FFTW], currently http://www.fftw.org/fftw-3.3.4.tar.gz .
* Extract the archive and enter the newly created directory.
* Configure:
```text
 ./configure  --prefix={{< basedir >}} CC="{{< cc >}}" CFLAGS="{{< cflags >}}" F77="{{< f90 >}}" F77FLAGS="{{< fflags >}}"
```
* Compile and install:
```text
 make
 make install
```

####  LibXC  

* Get {{< libxc >}}, currently http://www.tddft.org/programs/octopus/down.php?file=libxc/libxc-{{< libxc_version >}}.tar.gz .
* Extract the archive and enter the newly created directory.
```text
 tar -xvzf libxc-{{< libxc_version >}}.tar.gz
 cd libxc-{{< libxc_version >}}
```
* Configure:
```text
 ./configure --prefix={{< basedir >}} CC="{{< cc >}}" CFLAGS="{{< cflags >}}" FC="{{< f90 >}}" FCFLAGS="{{< fflags >}}"
```
* Compile and install:
```text
 make
 make install
```

###  Compilation of Octopus  

After compiling the libraries, now we are ready to compile {{< octopus >}}. These are the steps you have to follow:

* Download the last version of {{< octopus >}}: http://www.tddft.org/programs/octopus/down.php?file={{< octopus_version >}}/octopus-{{< octopus_version >}}.tar.gz
* Extract the file and enter the newly created directory.
* Define the following environment variables to be used by the configure script (we assume that you are using {{< name "bash" >}}, if you are using another shell, the commands should be modified accordingly):
```text
 export LIBS_BLAS={{< basedir >}}/lib/libblas.a
 export LIBS_LAPACK={{< basedir >}}/lib/liblapack.a
 export LIBS_FFT={{< basedir >}}/lib/libfftw3.a
```
* Now call the configure script:
```text
 ./configure CC="{{< cc >}}" CFLAGS="{{< cflags >}}" FC="{{< f90 >}}" FCFLAGS="{{< fflags >}}" --prefix={{< basedir >}} --with-gsl-prefix={{< basedir >}} --with-libxc-prefix={{< basedir >}}
```
* Compile and install
```text
 make
 make install
```
* If everything went fine, you should have {{< octopus >}} installed in the {{< basedir >}} directory, this means that the executables are in {{< basedir >}}{{< file "/bin/" >}}. You may want to add this last directory to your path to run octopus commands directly, otherwise you will have to give the full path.
* To test the compilation is correct, you can run the testsuite of {{< octopus >}}, that compares the results obtained with the created version against reference values. To do it, run
```text
 make check
```
All tests should be passed if the compilation is correct. If all of them fail, there is probably a problem with your executable, typically missing dynamic libraries. If
just some fail, there might be a more serious problem, we recommend you to look for help in the [[Mailing lists]].


{{< manual_foot prev="Manual:Updating to a new version" next="Tutorial:Running Octopus on Graphical Processing Units (GPUs)" >}}
---------------------------------------------
