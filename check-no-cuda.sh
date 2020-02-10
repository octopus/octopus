#!/bin/bash -l

module use ~/eb-buildbot/EasyBuild/modules/all

module purge

module load GSL/2.4-GCCcore-6.4.0
module load libxc/4.2.3-fosscuda-2018a

module load Octopus/bb-fosscuda-2018a
module load BerkeleyGW/1.2.0-fosscuda-2018a
module load NFFT/3.2.4-fosscuda-2018a
module load ETSF_IO/1.0.4-fosscuda-2018a
module load SPARSKIT2/2017.12.20-fosscuda-2018a
module load ELPA/2017.11.001-fosscuda-2018a
module load NLopt/2.4.2-fosscuda-2018a
module load libgd/2.2.5-GCCcore-6.4.0
module load libpspio/0.0.0-fosscuda-2018a
module load poke/ahi-fosscuda-2018a
module load libvdwxc/0.4.0-fosscuda-2018a



export INSTALLDIR=$HOME/Octopus_no_CUDA

export CFLAGS="$CFLAGS -g -O2 -march=native "
export FCFLAGS="-g -O2 -march=native -fbacktrace -fbounds-check -fcheck=all"
export CXXFLAGS="-O2 -march=native "

export CC=gcc
export FC=gfortran

export LIBS_BLAS="-L$EBROOTOPENBLAS/lib -lopenblas -pthread"

export CONFIG_FLAGS="\
    --with-libxc-prefix=$EBROOTLIBXC \
    --with-gsl-prefix=$EBROOTGSL \
    --with-poke-prefix=$EBROOTPOKE \
    --with-fftw-prefix=$EBROOTFFTW \
    --with-netcdf-prefix=$EBROOTNETCDFMINFORTRAN \
    --with-etsf-io-prefix=$EBROOTETSF_IO \
    --with-sparskit=$EBROOTSPARSKIT2/lib/libskit.a \
    --with-nlopt-prefix=$EBROOTNLOPT \
    --with-berkeleygw-prefix=$EBROOTBERKELEYGW \
    --with-pspio-prefix=$EBROOTLIBPSPIO \
    --with-libvdwxc-prefix=$EBROOTLIBVDWXC \
    --with-nfft=$EBROOTNFFT \
    --with-arpack=$EBROOTARPACKMINNG \
    --with-elpa-prefix=$EBROOTELPA \
    --enable-openmp \
    "

 
export OCT_TEST_NJOBS=8
export OMP_NUM_THREADS=1

make check > check.out 2> check.err
