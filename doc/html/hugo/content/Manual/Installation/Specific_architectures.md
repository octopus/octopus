---
title: "Specific architectures"
series: "Manual"
---


Here the [[Manual:Installation | general instructions]] for building Octopus are supplemented by instructions for some specific large supercomputers. Included also is how to configure the optional ETSF_IO library, and how you can run the testsuite in parallel.

#  Generic Ubuntu 10.10  
Install the following packages: 
```text
 * Required: gfortran, gcc, make, automake, m4, libtool, libgsl0-dev, libblas-dev, liblapack-dev, libfftw3-dev
 * Optional: mpi-default-dev, libgd2-xpm-dev, libsparskit-dev, libnetcdf-dev, libtrilinos-dev
```

[Unfortunately you should not use <tt>libblacs-mpi-dev, libscalapack-mpi-dev</tt> packages because they can give incorrect results. To use BLACS and ScaLAPACK (optional) you must build them yourself.]

First [http://www.tddft.org/programs/libxc/download download libxc], extract with {{< code "tar xzf" >}}, enter the resulting directory, and execute:

```text
 ./configure --prefix=`pwd` CC=gcc FC=gfortran FCCPP="/lib/cpp -ansi -freestanding" FCFLAGS="-O3 -ffree-line-length-none" CFLAGS=-O3
 make install
 make check
```

For octopus, download and extract, and execute the following (inserting the path where you put libxc):

```text
 ./configure --prefix=`pwd` CC=gcc CXX=g++ FC=gfortran FCCPP="/lib/cpp -ansi -freestanding" FCFLAGS="-O3 -ffree-line-length-none" CFLAGS=-O3 --with-libxc-prefix=''LIBXC_PATH''
 make install
```

If you are using MPI, replace the compilers by <tt>CC=/usr/bin/mpicc CXX=/usr/bin/mpicxx FC=/usr/bin/mpif90</tt> and add <tt>--enable-mpi</tt>.

To build ETSF_IO (optional), <tt>libnetcdf-dev</tt> is required.

```text
 ./configure --prefix=`pwd` CC=gcc FC=gfortran FCFLAGS="-O3 -ffree-line-length-none" CFLAGS="-O3" --with-netcdf-ldflags="-lnetcdff"
 make install
```

Add <tt>--with-etsf-io-prefix="$DIR/etsf_io-1.0.4"</tt> to the Octopus configure line, where <tt>$DIR/etsf_io-1.0.4</tt> is where you installed ETSF_IO.

(Aug 2018)

#  Generic MacOS  

Octopus is now available in MacPorts.
For more info, and how to build by hand with MacPorts libraries, see [[Compiling_octopus_in_OS_X]].

#  United States  

##  Sequoia/Vulcan (IBM Blue Gene/Q)  

Supercomputers at Lawrence Livermore National Laboratory based on the IBM Blue Gene/Q architecture. This was tested in Vulcan, but it should work on Sequoia as well.

https://computing.llnl.gov/tutorials/bgq/

```text
 export LIBS_BLAS="-L/usr/local/tools/essl/5.1/lib/ -lesslsmpbg"
 export LIBS_LAPACK="/usr/local/tools/lapack/lib/liblapack.a"
 export LIBS_FFT="-L/usr/local/tools/fftw-3.3.3/lib/ -lfftw3_omp"
 export FC_INTEGER_SIZE=4
 export CC_FORTRAN_INT=int
 export CC=mpixlc_r
 export FC=mpixlf95_r
 export LDFLAGS="-qsmp=omp"
 export CFLAGS="-g -O3"
 export FCFLAGS=$CFLAGS" -qxlf90=autodealloc -qessl"
 ./configure --with-libxc-prefix=$HOME --prefix=$HOME --host=powerpc32-unknown-linux-gnu  --build=powerpc64-unknown-linux-gnu --with-gsl-prefix=$HOME --enable-openmp --enable-mpi
```

(This build does not use Scalapack)

##  Stampede  

10 PFLOPS (PF) Dell Linux Cluster based on 6,400+ Dell PowerEdge server nodes, each outfitted with 2 Intel Xeon E5 (Sandy Bridge) processors and an Intel Xeon Phi Coprocessor (MIC Architecture). Texas Advanced Computing Center, University of Texas, Austin, USA.

https://portal.tacc.utexas.edu/user-guides/stampede

```text
 module swap mvapich2 impi
 module load fftw3 gsl hdf5 netcdf
 MKL_DIR=$MKLROOT
 ./configure --prefix=`pwd` CC=mpicc FC=mpif90 CFLAGS="-O3" FCFLAGS="-O3" --enable-mpi \
 --with-blas="-L$MKL_DIR/lib/intel64 -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group -lpthread" --with-fftw-prefix="$TACC_FFTW3_DIR" \
 --with-netcdf-prefix="$TACC_NETCDF_DIR" --enable-newuoa --with-libxc-prefix=`pwd`/../libxc-2.2.2 --with-blacs="$MKL_DIR/lib/intel64/libmkl_blacs_intelmpi_lp64.a" \
 --with-scalapack="$MKL_DIR/lib/intel64/libmkl_scalapack_lp64.a"
 make -j 6 install
```

For libxc:
```text
 ./configure --prefix=`pwd` CC=mpicc FC=mpif90 FCFLAGS="-O3" CFLAGS="-O3"
 make -j 6 install
 make check
```

testsuite script:

```text
 -!/usr/bin/env bash
 
 -SBATCH -n 16
 -SBATCH -p development
 -SBATCH -t 01:00:00
 -SBATCH --export=ALL
 
 module load perl
 WORKDIR=$PWD
 export TEMPDIRPATH=$SCRATCH/tmp
 cd $HOME/octopus
 export OCT_TEST_NJOBS=8
 export MPIEXEC=`which ibrun`
 make check &> $WORKDIR/makecheck
```

(13 September 2016)

##  Edison  
Cray XC30 at National Energy Research Scientific Computing Center (NERSC), Lawrence Berkeley National Laboratory, USA

http://www.nersc.gov/nusers/systems/

```text
 module load cray-fftw gsl cray-netcdf-hdf5parallel parpack
 
 ./configure --prefix=`pwd` --with-libxc-prefix=/usr/common/usg/libxc/3.0.0/intel/ivybridge --enable-mpi \
 CC=cc CXX=CC FC=ftn FCFLAGS="-fast -no-ipo" CFLAGS="-fast -no-ipo"  \
 FCCPP="/lib/cpp -ffreestanding" --with-fftw-prefix=$FFTW_DIR/..
 --with-arpack="$ARPACK" --with-parpack=no \
 --with-berkeleygw-prefix=/usr/common/software/berkeleygw/2.0
 make -j 6 install
```

To run the testsuite in parallel:

```text

-!/bin/bash -l
-SBATCH -J pulpo
-SBATCH -N 2
-SBATCH -p regular
-SBATCH -t 04:00:00
-SBATCH --export=ALL

cd $HOME/edison/octopus
export OMP_NUM_THREADS=1
export MPIEXEC=`which srun`
export TEMPDIRPATH=$SCRATCH/tmp
export OCT_TEST_NJOBS=15
make check &> $SLURM_SUBMIT_DIR/makecheck
```
</pre>

To build ''ETSF_IO'' 1.0.4 (optional): (some tests fail)

```text
 module load cray-netcdf-hdf5parallel
 ./configure --prefix=`pwd` FC=ftn FCFLAGS="-fast -no-ipo" CC=cc CFLAGS="-fast -no-ipo" \
 --with-netcdf-module-path="$NETCDF_DIR/include" --with-netcdf-prefix="$NETCDF_DIR"
 make -j 6 install
 make check
```

(version 8.1, July 2018)

##  Cori  
Cray XC40 at National Energy Research Scientific Computing Center (NERSC), Lawrence Berkeley National Laboratory, USA

http://www.nersc.gov/users/computational-systems/cori/

```text
 module load cray-fftw gsl cray-netcdf-hdf5parallel parpack
```

For libxc (you can use the installed version instead)
```text
 ./configure --prefix=`pwd` CC=cc FC=ftn FCFLAGS="-fast -no-ipo" CFLAGS="-fast -no-ipo"
 make -j 6 install
 make check
```

For octopus (note: using parpack here causes scalapack problems):
```text
 ./configure --prefix=`pwd` --with-libxc-prefix=/usr/common/software/libxc/4.2.3/intel/haswell --enable-mpi \
 CC=cc CXX=CC FC=ftn FCFLAGS="-fast -no-ipo" CFLAGS="-fast -no-ipo"  \
 FCCPP="/lib/cpp -ffreestanding" --with-fftw-prefix=$FFTW_DIR/.. \
 --with-arpack="$ARPACK" --with-parpack=no --with-berkeleygw-prefix=/usr/common/software/berkeleygw/2.0
 make -j 6 install
```

To run the testsuite in parallel:

```text

-!/bin/bash -l
-SBATCH -J pulpo
-SBATCH -n 64
-SBATCH -C haswell
-SBATCH -p regular
-SBATCH -t 04:00:00
-SBATCH --export=ALL

cd $HOME/cori/octopus
export OMP_NUM_THREADS=1
export MPIEXEC=`which srun`
export TEMPDIRPATH=$SCRATCH/tmp
export OCT_TEST_NJOBS=15
make check &> $SLURM_SUBMIT_DIR/makecheck
```
</pre>

To build ''ETSF_IO'' 1.0.4 (optional): (some tests fail)

```text
 module load cray-netcdf-hdf5parallel
 ./configure --prefix=`pwd` FC=ftn FCFLAGS="-fast -no-ipo" CC=cc CFLAGS="-fast -no-ipo"
 make -j 6 install
 make check
```

(30 July 2018)

##  Lonestar  
Dell Linux cluster at Texas Advanced Computing Center (TACC), University of Texas, Austin, USA

Part of National Science Foundation's [http://teragrid.org/ TeraGrid]

http://services.tacc.utexas.edu/index.php/lonestar-user-guide

```text
 module load gsl
 module load mkl
 module load netcdf
 module load fftw3
 cd libxc
 ./configure --prefix=`pwd`/.. --enable-mpi CC=mpicc FC=mpif90 F77=mpif90 FCFLAGS=-O3 CFLAGS=-O3
 make install
 cd ..
 export SLKPATH=/opt/apps/intel11_1/mvapich2_1_6/scalapack/1.8.0/lib
 ./configure --prefix=`pwd` --enable-mpi CC=mpicc FC=mpif90 \
 --with-blas="-Wl,--start-group $TACC_MKL_LIB/libmkl_intel_lp64.a $TACC_MKL_LIB/libmkl_sequential.a $TACC_MKL_LIB/libmkl_core.a -Wl,--end-group -lpthread" \
 --with-fft-lib="-L$TACC_FFTW3_LIB -lfftw3" --with-netcdf-prefix="$TACC_NETCDF_DIR" --disable-gdlib \
 --with-etsf-io-prefix="$HOME/etsf_io-1.0.3" --enable-newuoa --with-libxc-prefix=`pwd` \
 --with-blacs="$SLKPATH/libscalapack.a $SLKPATH/blacs_MPI-LINUX-0.a $SLKPATH/blacsF77init_MPI-LINUX-0.a $SLKPATH/blacs_MPI-LINUX-0.a"
 make install
```

To build ''ETSF_IO'':

```text
 module load netcdf
 ./configure --prefix=`pwd` FC=mpif90 --with-netcdf-module-path=$TACC_NETCDF_INC --with-netcdf-ldflags=-L$TACC_NETCDF_LIB FCFLAGS=-O3
 make install
 make check
```

To run the testsuite in parallel:

```text
 -!/bin/bash
 
 -$ -N test_pulpo
 -$ -cwd
 -$ -pe 6way 12
 -$ -q development
 -$ -l h_rt=00:45:00
 -$ -V
 
 cd $HOME/octopus
 export TEMPDIRPATH=$SCRATCH/tmp
 if [ ! -d $TEMPDIRPATH ]; then
     mkdir $TEMPDIRPATH
 fi
 
 export MPIEXEC=`which ibrun`
 export OCT_TEST_NPROCS=1
 make check &> makecheck
```

(2 May 2011)

##  Lawrencium  
Dell Linux cluster at Lawrence Berkeley National Laboratory, USA

Cluster available only to Lawrence Berkeley National Laboratory researchers

http://lrc.lbl.gov/html/Lawrencium.html

```text
 module load icc/10.1.018
 module load ifort/10.1.018
 module load mkl/2011.1.107
 module load openmpi/1.4.3-intel
 module load netcdf/4.0-intel
 module load fftw/3.1.2-intel
 module load autoconf/2.65
 module load gsl/1.13-intel
 module load gdlib/2.0.35
 module load libtool/2.2.6b
```

```text
 ./configure --prefix=`pwd` CC=mpicc FC=mpif90 CFLAGS="-O3" FCFLAGS="-O3 -check bounds" --enable-mpi \
 --with-blas="-L$MKLROOT/lib/intel64 -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group -lpthread" \
 --with-fft-lib="-L/global/software/centos-5.x86_64/modules/fftw/3.1.2-intel/lib -lfftw3" --with-netcdf-prefix="$NETCDF" --with-netcdf-include="$NETCDF/include" --enable-newuoa \
 --with-libxc-prefix=`pwd` --with-blacs="$MKL_DIR/lib/intel64/libmkl_blacs_openmpi_lp64.a" --with-scalapack="$MKL_DIR/lib/intel64/libmkl_scalapack_lp64.a"
```

To build ''ETSF_IO'':

```text
 module load netcdf/4.2-intel-p
 ./configure --prefix=`pwd` CC=mpicc FC=mpif90 CFLAGS="-O3" FCFLAGS="-O3" --with-netcdf-module-path=$NETCDF_DIR/include/ --with-netcdf-ldflags="-L$NETCDF_DIR/lib/ -lnetcdf"
```

(21 Aug 2012)

#  Europe  

##  Curie   

###  PFFT 1.0.5 and BigDFT 1.6.0  

The Curie supercomputer, owned by GENCI and operated into the TGCC by CEA, is the first French Tier0 system open to scientists through the French participation into the PRACE research infrastructure.

Curie is offering 3 different fractions of x86-64 computing resources for addressing a wide range of scientific challenges and offering an aggregate peak performance of 2 PetaFlops.

General information: http://www-hpc.cea.fr/en/complexe/tgcc-curie.htm

Compilation with PFFT version 1.0.5 and LibISF taken from BigDFT 1.6.0 (previously compiled):

```text
 -!/bin/bash
 
 module load c/intel
 module load fortran/intel
 module load bullxmpi
 module load mkl
 module load gsl/1.14
 
 WPREFIX=$WORKDIR/software
 PREFIX=$HOME/software
 export ISF_HOME="$WORKDIR/software/libisf"
 
 export CC=mpicc
 export CFLAGS="-xHost -O3 -m64 -ip -prec-div -static-intel -sox -I$ISF_HOME/include"
 export FC=mpif90
 export FCFLAGS="-xHost -O3 -m64 -ip -prec-div -static-intel -sox -I$PREFIX/fftw/3.3/include -I$ISF_HOME/include"
 
 export LIBS_BLAS="${MKL_LIBS}"
 export LIBS_FFT="-Wl,--start-group -L$PREFIX/pfft/1.0.5/lib -L$PREFIX/fftw/3.3/lib -lpfft -lfftw3 -lfftw3_mpi -Wl,--end-group"
 export LDFLAGS=$LDFLAGS" -L$ISF_HOME/lib -lisfsolver -lrt"
  
 ../configure --prefix=$WPREFIX/octopus/superciliosus_pfft \
   --disable-openmp --with-libxc-prefix=$WPREFIX/libxc/2.0.x \
   --disable-gdlib --with-gsl-prefix=/usr/local/gsl-1.14 \
   --enable-mpi --enable-newuoa \
   --with-pfft-prefix=$PREFIX/pfft/1.0.5
```

BigDFT was configured like this:

```text
 -!/bin/bash 
 
 export FC=mpif90
 
 ../configure --with-ext-linalg="-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lm" \
   --with-ext-linalg-path="-L/usr/local/Intel_compilers/c/composer_xe_2011_sp1.7.256/mkl/lib/intel64" \
   --prefix=$WORKDIR/software/bigdft
```


###  PFFT 1.0.7 and BigDFT 1.7.0  

```text
 -!/bin/bash
```

```text
 module load c/intel
 module load fortran/intel
 module load bullxmpi
 module load mkl
 module load gsl/1.14
 
 WPREFIX=$WORKDIR/poisson/software
 WLPREFIX=$WORKDIR/poisson/local
 export ISF_HOME="$WLPREFIX/libisf"
 export FFTW_HOME="$WLPREFIX/fftw-3.3.3"
 export PFFT_HOME="$WLPREFIX/pfft-1.0.7-alpha"
  
 export CC=mpicc
 export CFLAGS="-xHost -O3 -m64 -ip -prec-div -static-intel -sox -I$ISF_HOME/include"
 export FC=mpif90
 export FCFLAGS="-xHost -O3 -m64 -ip -prec-div -static-intel -sox -I$FFTW_HOME/include -I$ISF_HOME/include"
 
 export LIBS_BLAS="${MKL_LIBS}"
 export LIBS_FFT="-Wl,--start-group -L$PFFT_HOME/lib -L$FFTW_HOME/lib -lpfft -lfftw3 -lfftw3_mpi -Wl,--end-group"
 
 export LDFLAGS=" -L/usr/lib64 -L$ISF_HOME/lib -lPSolver-1 -lrt"
 export LD_LIBRARY_PATH=$ISF_HOME/lib:$LD_LIBRARY_PATH
 
 ../configure --prefix=$WPREFIX/octopus/superciliosus_pfft \
   --disable-openmp --with-libxc-prefix=$WPREFIX/libxc/11066 \
   --disable-gdlib --with-gsl-prefix=/usr/local/gsl-1.14 \
   --enable-mpi --enable-newuoa \
   --with-pfft-prefix=$PFFT_HOME
```

LibISF compilation:
```text
 -!/bin/bash
 
 module load c/intel
 module load fortran/intel
 module load bullxmpi
 module load mkl
 module load gsl/1.14
 
 export FC=mpif90
 
 ../configure --with-ext-linalg="-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lm" \
    --with-ext-linalg-path="-L/usr/local/Intel_compilers/c/composer_xe_2011_sp1.7.256/mkl/lib/intel64" \
    --prefix=$WORKDIR/poisson/local/libisf \
    --disable-libbigdft --disable-binaries --enable-dynamic-libraries
```

FFTW 3.3.3 compilation script is in http://www-user.tu-chemnitz.de/~mpip/software/install_fftw-3.3.3_gcc.sh and it is changed as follow:

```text
 < myprefix=$HOME/local
 --- 
 > module load c/intel
 > module load fortran/intel
 > module load bullxmpi
 > module load mkl
 > module load gsl/1.14
 > 
 > myprefix=$WORKDIR/poisson/local
 23c29,30
 < wget http://www.fftw.org/fftw-$FFTW_VERSION.tar.gz
 ---
 > - wget http://www.fftw.org/fftw-$FFTW_VERSION.tar.gz
 > cp ../fftw-$FFTW_VERSION.tar.gz .
 876c883
 < ./configure --prefix=$INSTALLDIR --enable-mpi --enable-openmp --enable-shared \
 ---
 > ./configure --prefix=$INSTALLDIR --enable-mpi --enable-openmp --enable-shared --disable-shared \ 
```


PFFTW 1.0.7 compilation script is in http://www-user.tu-chemnitz.de/~mpip/software/install_pfft-1.0.7-alpha_gcc.sh and it is changed as follow:

```text
 < myprefix=$HOME/local
 ---
 > 
 > module load c/intel
 > module load fortran/intel
 > module load bullxmpi
 > module load mkl
 > module load gsl/1.14
 > 
 > myprefix=$WORKDIR/poisson/local
 24c31,32
 < wget http://www.tu-chemnitz.de/~mpip/software/pfft-$PFFT_VERSION.tar.gz
 ---
 > -wget http://www.tu-chemnitz.de/~mpip/software/pfft-$PFFT_VERSION.tar.gz
 > cp ../pfft-$PFFT_VERSION.tar.gz .
```

##  Fermi/Juqueen  
IBM Blue Gene/Q, Cineca, Italy

General information: http://www.hpc.cineca.it/content/ibm-fermi-user-guide

The exactly the same script has been tested in the Juqueen, Forschungszentrum Jülich, Jülich Supercomputing Centre (JSC), Jülich, Germany. General information: http://www.fz-juelich.de/ias/jsc/EN/Expertise/Supercomputers/JUQUEEN/JUQUEEN_node.html

```text
 -/bin/bash
 
 module load gsl
 module load lapack
 module load blacs
 module load scalapack 
 
 export PFFT_HOME="$HOME/local/pfft-threads-1.0.5-alpha"
 export FFTW3_HOME="$HOME/local/fftw-threads-3.3.2"
 export LD_LIBRARY_PATH=$PFFT_HOME/lib:$FFTW3_HOME/lib:$LD_LIBRARY_PATH 
 
 export FC_INTEGER_SIZE=8
 export CC_FORTRAN_INT=int
 export CC=mpicc
 export CFLAGS="-g -Wl,-relax -O3 -I$PFFT_HOME/include -I$FFTW3_HOME/include"
 export FC=mpif90
 export FCFLAGS=$CFLAGS" -qxlf90=autodealloc -qessl -qsmp=omp "
 export LIBS_BLAS="-lesslsmpbg -L/opt/ibmmath/essl/4.4/lib -lesslbg"
 export LIBS_LAPACK="-L$LAPACK_LIB -llapack"
 
 export LIBS_FFT="-L$FFTW3_HOME/lib/ -lfftw3_mpi -lfftw3 -lm"
 echo "GSL DIR = $GSL_HOME "
 echo "PFFT DIR = $PFFT_HOME "
 export LDFLAGS=$LDFLAGS" -R/opt/ibmcmp/lib/bg/bglib \
    -L/bgsys/drivers/V1R4M2_200_2010-100508P/ppc/gnu-linux/powerpc-bgp-linux/lib \
    -L/bgsys/drivers/V1R4M2_200_2010-100508P/ppc/gnu-linux/lib/gcc/powerpc-bgp-linux/4.1.2 \
    -L/opt/ibmcmp/xlf/bg/11.1/bglib \
    -L/opt/ibmcmp/xlmass/bg/4.4/bglib \
    -L/opt/ibmcmp/xlsmp/bg/1.7/bglib \
    -L/bgsys/drivers/V1R4M2_200_2010-100508P/ppc/runtime/SPI \
    -L/bgsys/drivers/V1R4M2_200_2010-100508P/ppc/comm/sys/lib \
    -L/bgsys/drivers/V1R4M2_200_2010-100508P/ppc/comm/default/lib \
    -L/usr/local/zlib/v1.2.3/lib \
    -L/bgsys/drivers/ppcfloor/runtime/SPI \
    -L$PFFT_HOME/lib -L$FFTW3_HOME/lib \
    -lpfft -lfftw3_mpi -lfftw3 -lm \
    -ldl -lxl -lxlopt -lrt -lpthread"
 
 SVN_VERS=$(svn info .. | grep Revision | awk '{print $2}')
 echo $SVN_VERS
 
 ../configure --prefix=$HOME/software/octopus_omp \
    --with-libxc-prefix=$HOME/software/libxc_new \
    --with-fft-lib="$FFTW3_HOME/lib/libfftw3.a  $FFTW3_HOME/lib/libfftw3_omp.a" \
    --host=powerpc32-unknown-linux-gnu \
    --build=powerpc64-unknown-linux-gnu \
    --disable-gdlib \
    --disable-f90-forall \
    --with-blas="-L/opt/ibmmath/lib64 -lesslbg" \
    --with-lapack="-L$LAPACK_LIB -llapack" \
    --with-gsl-prefix=$GSL_HOME \
    --enable-mpi --enable-openmp
```

To build FFTW3. Small changes are made to the script taken from M. Pippig web site http://www-user.tu-chemnitz.de/~mpip/software.php?lang=en
It has a patch and is build with MPI and multithreaded (multithreaded may not work in Juqueen):

```text
 -!/bin/sh -e 
 
 myprefix=$HOME/local
 FFTW_VERSION="3.3.2"
 INSTALLDIR="$myprefix/fftw-threads-${FFTW_VERSION}"
 TMP="tmp-mine-fftw-${FFTW_VERSION}" 
 
 - bash check if directory exists
 if [ -d $TMP ]; then
         echo "Directory $TMP already exists. Delete it? (y/n)" 
 	read answer 
 	if [ ${answer} = "y" ]; then 
 		rm -rf $TMP 
 	else
 		echo "Program aborted."
 		exit 1
 	fi
  fi 
 
 mkdir $TMP && cd $TMP 
 
 wget http://www.fftw.org/fftw-$FFTW_VERSION.tar.gz
 gzip -dc fftw-$FFTW_VERSION.tar.gz | tar xvf -
 cd fftw-$FFTW_VERSION
 
 - fix bug in fftw alltoall routine that causes deadlocks
 patch -b -p0 <<\EOF
 --- mpi/transpose-alltoall.c
 +++ mpi/transpose-alltoall.c
 @@ -223,8 +223,8 @@
 	  sbo[pe] = (int) (pe * (b * p->tblock) * vn);
 	  rbs[pe] = (int) (db * bt * vn);
 	  rbo[pe] = (int) (pe * (p->block * bt) * vn);
 -	  if (sbs[pe] != (b * p->tblock) * vn
 -	      || rbs[pe] != (p->block * bt) * vn)
 +	  if (dbt != p->tblock
 +	      || db != p->block)
 	       equal_blocks = 0;
      }
      pln->send_block_sizes = sbs;
 EOF
 
 ./configure --build=powerpc64-bgq-linux-gnu --host=powerpc64-bgq-linux \
 --prefix=$INSTALLDIR --enable-mpi --enable-openmp --enable-threads --disable-shared \
 CC=mpixlc_r F77=mpixlf90_r MPICC=mpixlc_r MPICXX=mpixlcxx_r MPIFC=mpixlf90_r MPILIBS=" " \
 CFLAGS='-O3 -g' \
 FCFLAGS='-O3 -g' \ 
 
 make -j 4
 make install
```

To build PFFT, also taken from http://www-user.tu-chemnitz.de/~mpip/software.php?lang=en
```text
 -!/bin/sh -e 
 
 myprefix=$HOME/local
 PFFT_VERSION=1.0.5-alpha
 FFTW_VERSION=3.3.2
 INSTDIR=$myprefix/pfft-threads-$PFFT_VERSION
 FFTWDIR=$myprefix/fftw-threads-$FFTW_VERSION
 TMP="tmp-mine-pfft-$PFFT_VERSION" 
 
 - bash check if directory exists
 if [ -d $TMP ]; then
         echo "Directory $TMP already exists. Delete it? (y/n)"
 	read answer 
 	if [ ${answer} = "y" ]; then 
 		rm -rf $TMP 
 	else
 		echo "Program aborted."
 		exit 1
 	fi
 fi
 
 mkdir $TMP && cd $TMP
 
 wget http://www.tu-chemnitz.de/~mpip/software/pfft-$PFFT_VERSION.tar.gz
 gzip -dc pfft-$PFFT_VERSION.tar.gz | tar xvf -
 cd pfft-$PFFT_VERSION
 ./configure --build=powerpc64-bgq-linux-gnu --host=powerpc64-bgq-linux \
 --prefix=$INSTDIR --with-fftw3=$FFTWDIR \
 MPICC=mpif90 MPICXX=mpicxx MPIFC=mpif90 CC=mpicc FC=mpif90 \
 CFLAGS='-O3 -g -qmaxmem=-1' \
 FCFLAGS='-O3 -g -qmaxmem=-1' \
 --disable-shared
 
 make -j 4
 make install
```

##  Hydra   
Hydra is a supercomputer located at the RZG/MPG. More information: http://www.rzg.mpg.de/services/computing/hydra/about-the-system

Instruction that work with version 13307
```text
 
 -!/bin/bash
 
 . /etc/profile.d/modules.sh
 
 module purge  
 module load intel/14.0  
 module load mkl/11.1  
 module load mpi.ibm/1.3.0  
 module load gsl  
 module load fftw
 module load netcdf-serial
 
 export MKL="-L$MKL_HOME/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -lm"
 
 export NFF_HOME="/hydra/u/system/SLES11/soft/nfft/3.2.3/intel-14.0/mpi.ibm-1.3"
 export ISF_HOME=$HOME/local/PSolver/1.7.6
 
 export CC=mpiicc 
 export CFLAGS="-g -O3 -xHost -openmp -I$FFTW_HOME/include"
 export FC=mpiifort 
 export FCFLAGS="-g -O3 -xHost -openmp -I$FFTW_HOME/include" 
 
 export LDFLAGS="-I$FFTW_HOME/include -L$MKL_HOME/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -lpthread -lm -Xlinker -rpath=$MKL_HOME/lib/intel64:$GSL_HOME/lib:$NETCDF_HOME/lib:$NFF_HOME/lib" 
 
 export LD_LIBRARY_PATH=$ISF_HOME/lib:$LD_LIBRARY_PATH 
 
 ../configure --enable-mpi  --disable-gdlib --enable-openmp \
    --prefix=$HOME/software/octopus/repository \
    --with-gsl-prefix=$GSL_HOME \
    --with-nfft=$NFF_HOME  \
    --with-libxc-prefix=/hydra/u/system/SLES11/soft/libxc/2.0.3/intel-14.0/mpi.ibm-1.3 \
    --with-blas="$MKL" \
    --with-lapack="$MKL" \
    --with-isf-prefix=$ISF_HOME \
    --with-isf-include=$ISF_HOME/include \
    --with-metis-prefix=$HOME/local/parmetis \
    --with-parmetis-prefix=$HOME/local/parmetis \
    --with-netcdf-prefix=$NETCDF_HOME
```

##  Jugene  
FZJ-JSC IBM Blue Gene/P, JÜLICH SUPERCOMPUTING CENTRE (JSC), Germany

May also work on other Blue Gene/P computers

General information: http://www.fz-juelich.de/jsc/jugene

Usage information: http://www.fz-juelich.de/jsc/jugene/usage/
```text
 module load gsl;
 module load lapack;
 autoreconf -i;
 export FC_INTEGER_SIZE=4;
 export CC_FORTRAN_INT=int;
 export CC=mpixlc_r;
 export CFLAGS='-g -O3 -qarch=450d';
 export FC=mpixlf90_r;
 export FCFLAGS='-g -O3 -qarch=450d -qxlf90=autodealloc -qessl -qsmp=omp';
 export LIBS_BLAS="-lesslsmpbg -L/opt/ibmmath/essl/4.4/lib -lesslbg";
 export LIBS_LAPACK="-L$LAPACK_LIB -llapack";
 export LIBS_FFT="-L/bgsys/local/fftw3/lib/ -lfftw3 -lm";
 ./configure --prefix=$outdir \
      --host=powerpc32-unknown-linux-gnu \
      --build=powerpc64-unknown-linux-gnu \
      --disable-gdlib \
      --with-gsl-prefix=$GSL_DIR \
      --enable-mpi --enable-openmp;
```

To run you have to change to $WORK directory and run LoadLeveler. Below is an example of the tutorial (http://www.tddft.org/programs/octopus/wiki/index.php/Tutorial:Benzene_molecule) executing in 32 nodes in the SMP mode (4 threads per chip). The execution mode is hybrid; OpenMP/MPI. 

Here is a example job input script:
```text
 - @ job_name = Octopus_Sample_1
 - @ comment = "BGP Job by Size"
 - @ error = $(job_name).$(jobid).out
 - @ output = $(job_name).$(jobid).out
 - @ environment = COPY_ALL;
 - @ wall_clock_limit = 00:30:00
 - @ notification = error
 - @ notify_user = foo@gmail.com
 - @ job_type = bluegene
 - @ bg_size = 32
 - @ queue
 mpirun  -exe $OCTOPUS_HOME/bin/octopus_mpi -mode SMP -verbose 1 
```
To execute the above script:
```text
 llsubmit executing_script
```

##  MareNostrum II  
MareNostrum is a supercomputer in Europe, the most powerful in Spain. This is one of the seven supercomputers of the Spanish Supercomputing Network. The supercomputer consists of 2,560 JS21 blade computing nodes, each with 2 dual-core IBM 64-bit PowerPC 970MP processors running at 2.3 GHz for 10,240 CPUs in total. The computing nodes of MareNostrum communicate primarily through Myrinet.

You may need to compile GSL and FFTW libraries before starting then installation.

Here is a script to build Octopus:
```text
 -!/bin/bash
 
 export LD_LIBRARY_PATH=/usr/lib
 export CC=mpicc
 export FC=mpif90
 
 export CFLAGS="-q64 -O3 -qtune=ppc970 -qarch=ppc970 -qcache=auto -qnostrict -qignerrno -qinlglue"
 export FCFLAGS=$CFLAGS" -qessl -qxlf90=autodealloc" 
 
 export LIBS_BLAS="-L/usr/lib64 -lessl" 
 
 ./configure   \
    --with-lapack=/gpfs/apps/LAPACK/lib64/liblapack.a \
    --disable-gsltest \
    --disable-gdlib \
    --with-gsl-prefix=$outdirGSL \
    --with-fft-lib=$outdirFFT/lib/libfftw3.a  \
    --prefix=$outdir --enable-mpi --disable-f90-forall
```

As you can see in the script the --with-gsl-prefix has to point to your GSL and --with-fft-lib to your FFTW. You should compile those with the same CFLAGS and FCFLAGS. Here is an example of configuring GSL (we use GSL-1.11):
```text
 ./configure --prefix=$outdirGSL CC=gcc -m64 --enable-static --disable-shared
```
We use FFTW-3.1.3. To configure FFTW3 just write:
```text
 ./configure --prefix=$outdirFFT
```

Recently Marenostrum has changed to the SLURM manager. To execute a job you have to submit it with jobsubmit (more details in the [http://www.bsc.es/media/859.pdf User's Guide]).

##  MareNostrum III  

MareNostrum III is a supercomputer based on Intel SandyBridge processors, iDataPlex Compute Racks, a Linux Operating System and an Infiniband interconnection. More information: http://www.bsc.es/marenostrum-support-services/mn3

```text
 module load MKL
 module load GSL
 
 export PREFIX=$HOME/Software
 export FFTW3_HOME=$HOME/local/fftw-3.3.2
 export PFFT_HOME=$HOME/local/pfft-1.0.5-alpha
 
 export MKL_LIBS=$MKL_HOME/lib/intel64
 export MKLROOT=$MKL_HOME
 export MKL_LIBS="-L$MKLROOT/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 -lpthread -lm"
 export CC=mpicc
 export CFLAGS="-xHost -O3 -m64 -ip -prec-div -static-intel -sox -I$ISF_HOME/include  -openmp -I$MKLROOT/include"
 export FC=mpif90
 export FCFLAGS="-xHost -O3 -m64 -ip -prec-div -static-intel -sox -I$FFTW3_HOME/include  -openmp  -I$ISF_HOME/include -I$MKLROOT/include"
 
 export LIBS_BLAS="${MKL_LIBS}"
 export LIBS_FFT="-Wl,--start-group -L$PFFT_HOME/lib -L$FFTW3_HOME/lib -lpfft -lfftw3 -lfftw3_mpi -Wl,--end-group"
 
 
 ../configure --prefix=$PREFIX/octopus/ \
  --enable-openmp --with-libxc-prefix=$PREFIX/libxc \
  --disable-gdlib --with-gsl-prefix=/apps/GSL/1.15 \
  --enable-mpi --enable-newuoa \
  --with-pfft-prefix=$PFFT_HOME
```

##  Corvo  

Corvo is the computational cluster of the Nano-bio Spectroscopy Group. It consists of 1504 processor cores and 5284 GiB of RAM connected by an Infiniband network.

Script to build LIBFM library included in the Scafacos library (modified Scafacos svn version 1920, this is valid since version 9887 of Octopus):

```text
 ./configure --enable-fcs-solvers=fmm,direct --enable-fcs-fmm-comm=armci \
   --enable-fcs-fmm-unrolled --enable-fcs-fmm-max-mpol=40 \
   CFLAGS=-O3 CXXFLAGS=-O3 FCFLAGS=-O3 --disable-doc CXX=mpic++ CC=mpicc \
   FC=mpif90 --prefix=$HOME/Software/scafacos \
   LDFLAGS=-L/opt/scalapack/1.8.0/lib -L/opt/blacs/1.1/lib \
   -L/opt/gotoblas2/1.08/lib -L/opt/netcdf/4.0.1/lib -L/opt/lapack/3.2/lib \
   -L/opt/etsf_io/1.0.2/lib --enable-fcs-int=int --enable-fcs-float=double \
   --enable-fcs-integer=integer --enable-fcs-real=real*8 
 
 make
 make install
```

The scripts to compile FFTW 3.3.2 and PFFT 1.0.5

```text
 -!/bin/sh -e 
 
 myprefix=$HOME/local
 FFTW_VERSION="3.3.2"
 INSTALLDIR="$myprefix/fftw/${FFTW_VERSION}"
 TMP="tmp-fftw-${FFTW_VERSION}"
 
 - bash check if directory exists
 if [ -d $TMP ]; then
        echo "Directory $TMP already exists. Delete it? (y/n)"
        read answer
        if [ ${answer} = "y" ]; then
                rm -rf $TMP
        else
                echo "Program aborted."
                exit 1
        fi
 fi
 
 mkdir $TMP && cd $TMP
 
 wget http://www.fftw.org/fftw-$FFTW_VERSION.tar.gz
 gzip -dc fftw-$FFTW_VERSION.tar.gz | tar xvf -
 cd fftw-$FFTW_VERSION
 
 - fix bug in fftw alltoall routine that causes deadlocks
 patch -b -p0 <<\EOF
 --- mpi/transpose-alltoall.c
 +++ mpi/transpose-alltoall.c
 @@ -223,8 +223,8 @@
           sbo[pe] = (int) (pe * (b * p->tblock) * vn);
           rbs[pe] = (int) (db * bt * vn);
           rbo[pe] = (int) (pe * (p->block * bt) * vn);
 -         if (sbs[pe] != (b * p->tblock) * vn
 -             || rbs[pe] != (p->block * bt) * vn)
 +         if (dbt != p->tblock
 +             || db != p->block)
               equal_blocks = 0;
      }
      pln->send_block_sizes = sbs;
 EOF
 
 ./configure --prefix=$INSTALLDIR --enable-mpi --enable-openmp --enable-threads --disable-shared \
 CC="mpicc" F77="mpif90" MPICC="mpicc" \
 CFLAGS="-O3" FFLAGS="-O3" \
 MPILIBS=" "
 
 make -j 4
 make install
```


```text
 -!/bin/sh -e
 
 myprefix=/home/local
 PFFT_VERSION=1.0.5-alpha
 FFTW_VERSION=3.3.2
 INSTDIR=$myprefix/pfft/$PFFT_VERSION
 FFTWDIR=$myprefix/fftw/$FFTW_VERSION
 TMP="tmp-pfft-$PFFT_VERSION"
 module unload fftw/3.2.2
 
 - bash check if directory exists
 if [ -d $TMP ]; then
        echo "Directory $TMP already exists. Delete it? (y/n)"
        read answer
        if [ ${answer} = "y" ]; then
                rm -rf $TMP
        else
                echo "Program aborted."
                exit 1
        fi
 fi
 
 mkdir $TMP && cd $TMP
 export CC=mpicc
 export FC=mpif90 
 
 wget http://www.tu-chemnitz.de/~mpip/software/pfft-$PFFT_VERSION.tar.gz
 gzip -dc pfft-$PFFT_VERSION.tar.gz | tar xvf -
 cd pfft-$PFFT_VERSION
 ./configure --prefix=$INSTDIR --with-fftw3=$FFTWDIR --disable-shared
 
 make install
```

This is the scrip to compile Octopus with PFFT and FMM support (through Scafacos library):

```text
 module add gcc
 module add ifort
 module add gsl
 module add etsf_io
 module add lapack
 module add netcdf
 module add gotoblas2
 module add mpibull2 
 module add blacs
 module add scalapack
 
 export CC=mpicc
 export FC="mpif90"
 
 - It is important that "/opt/intel/Compiler//11.1/038/bin/intel64/" not to be in the PATH
 export PATH=/opt/mpi/mpibull2-1.3.9-14.s/bin:/opt/netcdf/4.0.1/bin:/opt/etsf_io/1.0.2/bin:/opt/gsl//1.13/bin:/opt/intel/composer_xe_2011_sp1.10.319//bin/intel64/:/opt/gcc//4.4/bin:/opt/slurm/bin:/usr/kerberos/bin:/opt/cuda//bin:/usr/local/bin:/bin:/usr/bin
 
 export LIBFM_HOME="$HOME/Software/scafacos"
 export PFFT_HOME="$HOME/local/pfft-1.0.5-alpha"
 export FFTW3_HOME="$HOME/local/fftw-3.3.2"
 export LIBS_PFFT="-L$PFFT_HOME/lib -L$FFTW3_HOME/lib -lpfft -lfftw3_mpi -lfftw3"
 export CFLAGS=" -I$PFFT_HOME/include -I$FFTW3_HOME/include -I$LIBFM_HOME/include"
 export LDFLAGS=$LDFLAGS" -L$LIBFM_HOME/lib -L$PFFT_HOME/lib -L$FFTW3_HOME/lib -lpfft -lfftw3_mpi -lfftw3 -lm "
 export FCFLAGS=$CFLAGS
 export CPPFLAGS=" -I$LIBFM_HOME/include "
 export LD_LIBRARY_PATH=$LIBFM_HOME/lib:$LD_LIBRARY_PATH
  
 ../configure --enable-mpi --with-libxc-prefix=$HOME/Software/libxc_9882 \
    --prefix=$HOME/Software/octopus \
    --with-blacs=/opt/blacs/1.1/lib/libblacs.a --with-scalapack \
    --with-libfm="-L$LIBFM_HOME/lib -lfcs4fortran -lfcs -lfcs_direct -lfcs_fmm -lfcs_near -lfcs_gridsort -lfcs_common" \
    --disable-openmp
```
##  LaPalma  

LaPalma is a supercomputer of the IAC. It is a part of the old MareNostrum II computer. It was not possible to compile with XL compilers, so GNU 4.6.1 version is used.

Previously to Octopus, Libxc (2.2), FFTW (3.3.4), GSL (1.16) and OpenBLAS (0.2.13) were required to compile.

GSL
```text
 export LD_LIBRARY_PATH=/usr/lib:/gpfs/apps/MPICH2/slurm-2.5.6/64/lib/:/gpfs/apps/GCC/4.9.2/lib64/
 export CC="/gpfs/apps/MPICH2/mx/default/64/bin/mpicc -cc=/gpfs/apps/GCC/4.9.2/bin/gcc"
 export FC="/gpfs/apps/MPICH2/mx/default/64/bin/mpif90 -fc=/gpfs/apps/GCC/4.9.2/bin/gfortran"
 
 ./configure --prefix=$HOME/local/gsl/gnu --enable-static --disable-shared
```

OpenBLAS
```text
 make PREFIX=$HOME/local/openblas/4.6.1 BINARY=64 CC=/gpfs/apps/GCC/4.6.1/bin/gcc FC=/gpfs/apps/GCC/4.6.1/bin/gfortran
```
Libxc
```text
 export PATH=/gpfs/apps/MPICH2/mx/default/64/bin:$PATH
 export LD_LIBRARY_PATH=/gpfs/apps/MPICH2/mx/default/64/lib:/gpfs/apps/MPICH2/slurm/64/lib
 
 export LD_LIBRARY_PATH=/gpfs/apps/GCC/4.6.1/lib64/:$HOME/local/openblas/lib:$LD_LIBRARY_PATH
 export CC="/gpfs/apps/MPICH2/mx/default/64/bin/mpicc -cc=/gpfs/apps/GCC/4.6.1/bin/gcc"
 export FC="/gpfs/apps/MPICH2/mx/default/64/bin/mpif90 -fc=/gpfs/apps/GCC/4.6.1/bin/gfortran"
 
 export CFLAGS="-m64 -O3"
 export FCFLAGS=$CFLAGS
 
 ./configure --prefix=$HOME/local/libxc/gnu/4.6.1/2.2
```

Octopus 

```text
 export PATH=/gpfs/apps/MPICH2/mx/default/64/bin:$PATH
 export LD_LIBRARY_PATH=/gpfs/apps/MPICH2/mx/default/64/lib:/gpfs/apps/MPICH2/slurm/64/lib
 
 export LD_LIBRARY_PATH=/gpfs/apps/GCC/4.6.1/lib64/:$HOME/local/openblas/4.6.1/lib:$LD_LIBRARY_PATH
 export CC="/gpfs/apps/MPICH2/mx/default/64/bin/mpicc -cc=/gpfs/apps/GCC/4.6.1/bin/gcc"
 export FC="/gpfs/apps/MPICH2/mx/default/64/bin/mpif90 -fc=/gpfs/apps/GCC/4.6.1/bin/gfortran"
 
 export CFLAGS="-m64 -O2 -I $PWD/external_libs/isf/wrappers "
 export FCFLAGS=$CFLAGS
 
 ../configure   \
    --with-blas=$HOME/local/openblas/4.6.1/lib/libopenblas.so.0 \
    --with-lapack=/gpfs/apps/LAPACK/lib64/liblapack.a \
    --disable-gdlib \
    --with-gsl-prefix=$HOME/local/gsl/gnu \
    --with-fft-lib=$HOME/local/fftw/3.3.4/lib/libfftw3.a  \
    --with-libxc-prefix=$HOME/local/libxc/gnu/4.6.1/2.2 \
    --prefix=/gpfs/projects/ehu31/software/octopus/gnu/4.6.1 \
    --enable-mpi --disable-f90-forall
```
##  XBES  

Located at the University of Hamburg in Hamburg. It has a login node and a 24 nodes with Intel Xeon. Here a basic configuration script:

```text
 ./configure \
 PATH="$PATH:/home/common/lib/gsl-1.16-intel/bin" \
 CC="mpiicc -m64" CFLAGS=" -O3 -march=native  -I${MKLROOT}/include/intel64/lp64 -I${MKLROOT}/include" \
 FC="mpiifort -m64" FCFLAGS="-O3  -I${MKLROOT}/include/intel64/lp64 -I${MKLROOT}/include" \
 GSL_CONFIG=/home/common/lib/gsl-1.16-intel/bin/gsl-config --with-gsl-prefix="/home/common/lib/gsl-1.16-intel/" \
 --with-libxc-prefix=$LIBXC_PREFIX   \
 --with-blas="${MKLROOT}/lib/intel64/libmkl_blas95_lp64.a ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm -ldl" \
 --with-lapack="${MKLROOT}/lib/intel64/libmkl_blas95_lp64.a ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm -ldl" \
 --with-fft-lib=/home/common/lib/fftw-3.3.4-intel/lib/libfftw3.a
```

##  Magerit (CeSViMa)  

Magerit is a cluster located in "Centro de Supercomputación y Visualización de Madrid".  It consists on 245 nodes eServer BladeCenter PS702 with 2 Power7 processors each of 8 cores at 3'3 GHz (422,4 GFlops) and with 32 GB de RAM.
Here the "optimal" configuration script which uses GNU compilers with ESSL IBM library.

```text
 module purge
 module load gcc/4.4
 module load openmpi/1.6.3
 
 export LIBSDIR='/sw/openmpi/'
 
 export SOFTWARE_DIR="$HOME/SOFTWARE/exec"
 
 export FC=mpif90
 export CC=mpicc
 export CXX=mpicxx
 
 export STDFLAGS="-O3 -mcpu=power7 -mtune=power7"
 
 export FCFLAGS=$STDFLAGS
 export CFLAGS=$STDFLAGS
 export CXXLAGS=$STDFLAGS
 
 export C_TYPE="gnu64-4.4-essl"
 
 export branch=`awk '{print $2}' ../.git/HEAD`
 export branch=`basename $branch`
 
 INSTALLDIR="$SOFTWARE_DIR/octopus/git-$C_TYPE/$branch"
 
 -GSL  1.16 --compiled with gcc 4.4.6
 export GSL_HOME="$LIBSDIR/GSL/1.16/"
 
 -LAPACK gcc 4.4
 export LAPACK_HOME="$LIBSDIR/LAPACK/3.4.2-gnu64-4.4"
 
 -PARMETIS gcc 4.7.2
 export PARMETIS_HOME="$LIBSDIR/METIS/PARMETIS-4.0.3/"
 
 -METIS gcc 4.7.2
 export METIS_HOME="$LIBSDIR/METIS/METIS-5.1.0/"
 
 -FFTW3.3.4 gcc 4.4
 export FFTW_HOME="$LIBSDIR/FFTW/3.3.4-gnu64-4.4"
 export LD_LIBRARY_PATH=$FFTW_HOME/lib:$LD_LIBRARY_PATH
 
 -SCALAPACK gcc  4.4.6
 export SCALAPACK_HOME="$LIBSDIR/SCALAPACK/2.0.2-gnu64-4.4/"
 export LD_LIBRARY_PATH=$SCALAPACK_HOME/lib:$LD_LIBRARY_PATH
 
 -LIBXC
 export LIBXC_HOME="$SOFTWARE_DIR/libxc/3.0.0-$C_TYPE/"
 
 -LIBVDWXC
 export LIBVDWXC_HOME="$SOFTWARE_DIR/libvdwxc/0.2.0-$C_TYPE/"
 export CPPFLAGS="-I$LIBVDWXC_HOME/include -DHAVE_LIBVDWXC"
 -export LIBS=-lvdwxcfort
 export LDFLAGS="-L$LIBVDWXC_HOME/lib -lvdwxcfort -Wl,-rpath=$LIBVDWXC_HOME/lib"
 export LD_LIBRARY_PATH=$LIBVDWXC_HOME/lib:$LD_LIBRARY_PATH
 
 export LD_LIBRARY_PATH=/opt/ibmcmp/lib64/:/opt/ibmcmp/xlsmp/2.1/lib64/:/opt/ibmcmp/xlf/13.1/lib64/:$LD_LIBRARY_PATH
 
 export LIBS_BLAS=" -L$SCALAPACK_HOME/lib -L$LAPACK_HOME/lib -lscalapack -llapack -lessl -lblacssmp -L$FFTW_HOME/lib -lfftw3 -lfftw3_threads -lfftw3_mpi / opt/gnu/gcc/4.4.6/lib64/libgfortran.a -L/opt/ibmcmp/xlf/13.1/lib64/ -lxlf90_r -lxlfmath -lxl -L/opt/ibmcmp/xlsmp/2.1/lib64/ -lxlomp_ser -Wl,--rpath /opt/ibmcmp/xlf/13.1/lib64/  /opt/ibmcmp/xlf/13.1/lib64/libxl.a -R/opt/ibmcmp/lib64/"
 
 make distclean
 make clean
 ../configure --prefix=$INSTALLDIR \
     --disable-gdlib --disable-gsltest \
     --with-libxc-prefix=$LIBXC_HOME \
     --with-fftw-prefix=$FFTW_HOME \
     --with-gsl-prefix=$GSL_HOME \
     --with-metis-prefix=$METIS_HOME\
     --with-parmetis-prefix=$PARMETIS_HOME \
     --enable-mpi \
     --enable-openmp
```


##  MPCDF systems (draco, cobra, raven)  

Draco, Cobra and Raven are the HPC machines of the Max-Planck Compute and Data Facility in Garching. On these, Octopus can be built with the following script:

```text
 -! /bin/bash -l
 - compilation script for octopus
 - needs to be executed in a subfolder (e.g. _build) of the root directory
 -
 - if you want to skip the configure step, call as ./build_octopus.sh noconfig
 - if you want to skip the configure step and the dependency checking,
 -   call as ./build_octopus.sh noconfig nodep
 
 if [[ ! -f ../configure.ac ]]; then
   echo "Error! Please execute this script in a subfolder of the root directory."
   exit 1
 fi
 
 compiler=intel
 -compiler=gnu
 cuda=yes
 
 echo Using $compiler compiler.
 [[ $cuda == yes ]] && echo Building with CUDA support.
 
 module purge
 if [[ $compiler == intel ]]; then
   - newest intel version
   module load intel/19.1.2 impi/2019.8 mkl/2020.2
   export CC=mpiicc
   export FC=mpiifort
   export CXX=mpiicpc
   export CFLAGS="-O3 -xCORE-AVX512 -qopt-zmm-usage=high -fma -ip"
   export FCFLAGS="$CFLAGS"
   export CXXFLAGS="$CFLAGS"
   export MKL="-L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 -liomp5 -lpthread -lm -ldl"
 elif [[ $compiler == gnu ]]; then
   module load gcc/9 impi/2019.8 mkl/2020.2
   export CC=mpicc
   export FC=mpif90
   export CXX=mpicxx
   export CFLAGS="-O3 -march=skylake-avx512 -g"
   export FCFLAGS="$CFLAGS"
   export CXXFLAGS="$CFLAGS"
   export MKL="-L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 -lgomp -lpthread -lm -ldl"
 else
   echo "Compiler $compiler unknown."
   exit 1
 fi
 
 module load gsl hdf5-serial netcdf-serial libxc metis parmetis
 export LDFLAGS="-Xlinker -rpath=$MKL_HOME/lib/intel64:$GSL_HOME/lib:$NETCDF_HOME/lib:$ELPA_HOME/lib:$METIS_HOME/lib:$PARMETIS_HOME/lib"
 
 if [[ $cuda == yes ]]; then
   module load cuda/10.2
   CUDA_FLAGS="--enable-cuda --enable-nvtx --with-cuda-prefix=$CUDA_HOME"
   LDFLAGS="$LDFLAGS:$CUDA_HOME/lib64"
 else
   CUDA_FLAGS=""
 fi
 module list
 
 export INSTALLDIR=$PWD/installed
 
 if [[ "$1" != "noconfig" ]]; then
   pushd .. && autoreconf -i && popd  
   ../configure $CUDA_FLAGS \
     FCFLAGS_FFTW="-I$MKLROOT/include/fftw" \
     FCCPP="cpp -ffreestanding" \
     --prefix=$INSTALLDIR \
     --enable-mpi --enable-openmp \
     --disable-gdlib \
     --with-gsl-prefix="$GSL_HOME" \
     --with-libxc-prefix="$LIBXC_HOME" \
     --with-blas="$MKL" \
     --with-lapack="$MKL" \
     --with-blacs="$MKL" \
     --with-scalapack="$MKL" \
     --with-netcdf-prefix="$NETCDF_HOME" \
     --with-metis-prefix="$METIS_HOME" \
     --with-parmetis-prefix="$PARMETIS_HOME" \
     || exit 1
 fi
 
 echo "\n\nBuilding octopus...\n"
 
 if [[ "$2" == "nodep" ]]; then
   make NODEP=1 -j20 && make NODEP=1 install || exit 1
 else
   make -j20 && make install || exit 1
 fi
 
 mkdir -p $INSTALLDIR/.build.doc/
 cp -f config.log $INSTALLDIR/.build.doc/
 cp -f $0 $INSTALLDIR/.build.doc/
 
 echo "... done"
```

{{< manual_foot prev="Tutorial:Running Octopus on Graphical Processing Units (GPUs)" next="Manual:Appendix:Porting Octopus and Platform Specific Instructions" >}}
---------------------------------------------
