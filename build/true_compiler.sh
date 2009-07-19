#!/bin/sh
# 
# $Id: true_compiler.sh 4655 2008-10-13 17:31:10Z acastro $

# If the argument is mpicc or mpif90, prints the true compiler that is
# being called. For the moment it only works if the -show argument is
# accepted by the wrapper (mpich, openmpi and derivatives do).

if [ x$1 == xmpicc -o x$1 == xmpif90 ]
then
echo -n "("`$1 -show | cut -f 1 -d" "`")"
fi
