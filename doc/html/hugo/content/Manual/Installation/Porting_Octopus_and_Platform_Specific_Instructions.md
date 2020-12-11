---
title: "Porting Octopus and Platform Specific Instructions"
series: "Manual"
---


This page contains information about {{< octopus >}} portability, with specific information to compile and run octopus for many architectures. If you managed to compile octopus for a different system, please contribute.
Warning: this information is quite out of date and may no longer be valid.

###  General information and tips about compilers  

* {{< octopus >}} is developed in the {{< name "GNU" >}} environment and sometimes we depend on system features that are {{< name "GNU" >}} extensions without knowing it. These are bugs and we will try to fix them; please report any problem that you find.

* If you have problems with the C compiler, try to use {{< name "gcc" >}}. It normally works better than vendor compilers and it's available on most systems. However, be careful with locally installed versions of gcc: sometimes they don't work.

* The {{< name "Fortran" >}} {{< code "//" >}} concatenation operator is sometimes recognized as a {{< name "C++" >}}-style comment and the preprocessor gives erroneous results: sometimes it doesn't expand macros after it or simply eliminates what comes after. To solve this problem, use the preprocessor with the {{< name "-C" >}} (keep comments) and {{< name "-ansi" >}} or equivalent options (in {{< name "ANSI C" >}} {{< code "//" >}} is not a comment).

* If you are compiling in dual 32/64-bit architectures like {{< name "PowerPC" >}}, {{< name "UltraSparc" >}} or {{< name "AMD64" >}} systems here are some tips: 
** A 64-bit version of {{< octopus >}} is only needed if you are going to use more than 2-3 {{< name "Gb" >}} of physical {{< name "RAM" >}}.
** Some operating systems have 64 bits {{< name "kernels" >}} and 32 bits {{< name "userspace" >}} ({{< name "Solaris" >}}, {{< name "OS X" >}}); if you want a 64-bit {{< octopus >}} there, you have to compile all required libraries in 64-bit (normally a 64-bit {{< name "libc" >}} is available).
** Typically Linux distributions {{< emph "for" >}} {{< name "AMD64" >}} have a 64-bit {{< name "userspace" >}}, so you will get a 64-bit executable there.

###  SSE2 support  

* We have some SSE2 code written using compiler primitives that can give an important increase in perfomance. For this you need hardware support (AMD Opteron/Athlon 64/Sempron/Turion or Intel Pentium 4 or newer) and compiler support, supported compilers are GCC and pathcc. For gcc you need to put the correct {{< flag "-march" >}} flags (for example {{< code "<nowiki>-march=opteron</nowiki>" >}} or {{< code "<nowiki>-march=pentium4</nowiki>" >}}).

* Besides this, for {{< name "x86" >}} (32 bits) you have to link dynamically because we have to use a tweaked malloc function that doesn't work with static linking. For x86_64 (64 bits) this is not needed.

###  Operating systems  

####  Linux  

The main development operating system for {{< octopus >}}.

####  Solaris  

Octopus compiles correctly either with sun compilers or gcc/gfortran. By default {{< name "Solaris" >}} doesn't have {{< name "GNU coreutils" >}}, so some test won't run.

####  Tru 64  

It works.

####  Mac OS X  

It works. Don't try to compile static binaries, they are not supported by the OS.

####  Windows  

Toy operating systems are not supported for the moment, sorry.

###  Compilers  

####  Intel Compiler for x86/x86_64  

* status: ok
* Version 9 and version 7.1 Build 20040901Z are ok. Versions 8 and 8.1 can be problematic.
* Recommended flags: FCFLAGS="-u -zero -fpp1 -nbs -pc80 -pad -align -unroll -O3 -ip -tpp7 -xW"
* {{< name "Intel" >}} artificially blocks their compilers from using certain optimization in non-{{< name "Intel" >}} processors.
* With Octopus 3.2.0, use of the flags <tt>-check all -traceback</tt> with <tt>ifort</tt> 10.1.018 will cause an internal compiler segmentation fault while compiling <tt>src/grid/mesh_init.F90</tt>.

####  Intel Compiler for Itanium  

* status: ok
* Version: 8.1.033 (older 8 releases and version 9 are reported to cause problems), version 10 works but it is much slower than 8.1.
* Recommended flags:
** FCFLAGS="-O3 -tpp2 -ip -IPF_fp_relaxed -ftz -align all -pad"
** CFLAGS="-O3 -tpp2 -ip -IPF_fp_relaxed -ftz"

####  [http://www.open64.net/ Open64]  

This is an open source compiler based on the liberated code of SGI MIPSpro compiler. It is available for x86, x86_64 and Itanium architectures.

####  [http://www.pathscale.com/ Pathscale Fortran Compiler] 

* Versions tested: 2.2, 2.3 and 2.5
* Architecure: x86, AMD64
* Recommended flags:
```text

FCFLAGS="-Wall -O3 -march=auto -mcpu=auto -OPT:Ofast -fno-math-errno"
```
</pre>
* Issues: 
** Everything works.
** It's necessary to compile blas/lapack with the same compiler.

####  NAG compiler  

AMD64:
```text

FCFLAGS="-colour -kind=byte -mismatch_all -abi=64 -ieee=full -O4 -Ounroll=4"
```
</pre>

####  [http://gcc.gnu.org/ GNU C Compiler (gcc)]  

####  [http://gcc.gnu.org/fortran/ GNU Fortran (gfortran)]  

* Status: ok.
* Version: gcc version 4.1.1 or newer. (4.1.0 {{< emph "does not" >}} work) For the parallel version you need at least gfortran 4.3.
* You may also need to compile blas, lapack and fftw3 using that specific gfortran version.
* Some recommended flags: <tt>-march=athlon64 -msse2 -mfpmath=sse -malign-double -funroll-loops -O3</tt>

#### [http://www.g95.org/ g95] 

* Status: works
* Tested architectures: x86/Linux, PowerPC/Darwin
* Version: version 4.0.3 (g95 0.91!) May 24 2007
* G95 doesn't recognize the linemarkers created by the preprocessor, so it's necessary to pass the -P flag to cpp.
* Flags:
```text

FC=g95
FCFLAGS="-O3 -funroll-loops -ffast-math"
FCCPP="cpp -ansi-P"
```
</pre>

There may be problems with versions 0.92 or 0.93, depending on the underlying version of gcc. See [[G95]] for info on building version 0.94 with gcc 4.2.4.

####  Portland 6  

Flags:

```text

FCFLAGS="-fast -mcmodel=medium -O4"
```
</pre>

Known problems:

The following problem with the PGI compiler version 6.0 and MPICH version 1.2.6 on {{< name "x86_64" >}} has been reported:

The MPI detection during the {{< command "configure" >}} step does not work properly. This may lead to compilation failures on e. g. the file {{< file "par_vec.F90" >}}. This problem is considered a bug in either the PGI compiler or the MPICH implementation. Please apply the following change by hand after running {{< command "configure" >}}:

In the file {{< file "config.h" >}}, replace the line 

```text

/* -undef MPI_H */
```
</pre>

by

```text

-define MPI_H 1
```
</pre>

and remove the line

```text

-define MPI_MOD 1
```
</pre>

####  Portland 7, 8, 9  

Flags (tested on Cray XT4):

```text

FCFLAGS="-O2 -Munroll=c:1 -Mnoframe -Mlre -Mscalarsse -Mcache_align -Mflushz"
```
</pre>

The configure script may fail in the part <tt>checking for Fortran libraries of mpif90</tt> for <tt>autoconf</tt> version 2.59 or earlier. The solution is to update <tt>autoconf</tt> to 2.60 or later, or manually set <tt>FCLIBS</tt> in the <tt>configure</tt> command line to remove a spurious apostrophe.

####  Portland 10  

For Octopus 3.2.0, the file <tt>src/basic/lookup.F90</tt> is incorrectly optimized yielding many segmentation faults in the testsuite. With PGI 10.5 the optimization flag should be <tt>-O2</tt> or less; with PGI 10.8 the optimization flag should be <tt>-O1</tt> or less. Note that <tt>-fast</tt> and <tt>-fastsse</tt> are between <tt>-O2</tt> and <tt>-O3</tt>. For later versions of Octopus, a PGI pragma compels this file to be <tt>-O0</tt> regardless of what is specified in <tt>FCFLAGS</tt>, so you may safely set <tt>FCFLAGS</tt> to <tt>-fast</tt>.

####  Portland 11  

11.4 does not work and will crash with glibc memory corruption errors. 11.7 is fine.

####  Portland 12  

12.5 and 12.6 cannot compile due to an internal compiler errors of this form:

```text
 PGF90-S-0000-Internal compiler error. sym_of_ast: unexpected ast    6034 (simul_box.F90: 1103)
```

12.4 and 12.9 are ok.

####  Absoft  
Flags x86:
```text

FCFLAGS="-O3 -YEXT_NAMES=LCS -YEXT_SFX=_"
```
</pre>

Flags amd64/em64t:
```text

FCFLAGS="-O3 -mcmodel=medium -m64 -cpu:host -YEXT_NAMES=LCS -YEXT_SFX=_"
```
</pre>

####  Compaq compiler 

```text

FCFLAGS="-align dcommons -fast -tune host -arch host -noautomatic"
```
</pre>

####  Xlf  
* Status: works
-bmaxdata:0x80000000 -qmaxmem=-1 -qsuffix=f=f90 -Q -O5 -qstrict -qtune=auto -qarch=auto -qhot -qipa

* Because of the exotic mixture of {{< name "MAC OS" >}} and {{< name "BSD" >}}, this system is not very standard. Compiling {{< octopus >}} can be problematic.

* {{< name "OS X" >}} doesn't support static linking of binaries, so don't try.

####  SGI MIPS 

-O3 -INLINE -n32 -LANG:recursive=on

####  Sun Studio  

You can download this compiler for free, it supports Linux and Solaris over x86, amd64 and sparc. A very fast compiler but quite buggy.

* Flags:
** CFLAGS="-fast -xprefetch -xvector=simd -D__SSE2__"
** FCFLAGS=$FLAGS

###  MPI Implementations  

####  OpenMPI  

####  MPICH2  

####  [http://www.sgi.com/products/software/mpt/ SGI MPT]  

####  Intel MPI  

####  [http://www.sun.com/software/products/clustertools/ Sun HPC ClusterTools] 

####  [http://mvapich.cse.ohio-state.edu/ MVAPICH]  

###  NetCDF  

{{< octopus >}} uses the Fortran 90 interface of {{< name "netCDF" >}}, this means that it's likely that you will have to compile it using the same compiler you will use to compile {{< octopus >}}. You can get the sources and follow installation instructions from the [http://www.unidata.ucar.edu/software/netcdf/ {{< name "NetCDF" >}} site].

###  BLAS and LAPACK  

These are standard libraries that provide a series of common vector and matrix operations. {{< octopus >}} uses as much as possible this libraries. There are several version available depending on your hardware. Around 40% of {{< octopus >}} execution time is spend in BLAS level 1 routines, so getting a fast implementation for your hardware might be important. On the other hand, Lapack performance is not very important.

####  AMD ACML  

This is the AMD Mathematical Library optimized to run in Athlon and Opteron processors. You can get a free copy from http://developer.amd.com/acml.jsp .

#### [http://math-atlas.sourceforge.net/ ATLAS]  

####  Compaq CXML  

#### [GOTO BLAS]  

Probably the fastest implementation of blas, source code is available and it can be compiled in many architectures.

####  Intel MKL  

See https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor for MKL's advice on the proper way to link. Here is an example, in which {{< code "--with-lapack" >}} is left blank because it is included in {{< code "--with-blas" >}}.

```text
 MKL_DIR=/opt/intel/mkl/lib/lintel64
 --with-blas="-L$MKL_DIR -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group -lpthread"
 --with-blacs="$MKL_DIR/libmkl_blacs_intelmpi_lp64.a" --with-scalapack="$MKL_DIR/libmkl_scalapack_lp64.a"
```

####  Netlib   

The reference implementation of BLAS and Lapack. It is available in most linux distributions. You can get the source code from http://www.netlib.org/blas http://www.netlib.org/lapack .

{{< manual_foot prev="Manual:Specific architectures" next="Manual:Appendix:Reference Manual" >}}
---------------------------------------------
