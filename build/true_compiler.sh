#!/usr/bin/env bash
# 
# $Id: true_compiler.sh 4655 2008-10-13 17:31:10Z acastro $

# If the argument is mpicc or mpif90, prints the true compiler that is
# being called. For the moment it only works if the -show argument is
# accepted by the wrapper (mpich, openmpi and derivatives do).
# It also works for the Cray compiler wrappers cc and ftn.

RESULT="(unknown)"

# MPI
if echo $1 | grep mpi > /dev/null; then
  RESULT="("`$1 -show | cut -f 1 -d" "`")"
# Cray compiler wrappers
elif [ x$1 == xcc -o x$1 == xftn -o x$1 == xCC ]; then
  if $1 -show &> /dev/null; then
      RESULT="("`$1 -show | grep DRIVERNAME | cut -f 2 -d"="`")"
  fi
fi

echo -n "$RESULT"
