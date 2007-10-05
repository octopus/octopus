project (LibNBC)

SET(LibNBC_SRC nbc.c nbc_op.c dict.c hb_tree.c nbc_fortran.c)
SET(LibNBC_COLLALGS_SRC nbc_ialltoallv.c nbc_ibarrier.c nbc_ibcast.c nbc_igather.c nbc_igatherv.c nbc_iscatter.c nbc_iscatterv.c nbc_iallgather.c nbc_iallgatherv.c nbc_ireduce.c nbc_ireduce_scatter.c nbc_iallreduce.c nbc_iscan.c nbc_ialltoall.c)

SET(CMAKE_VERBOSE_MAKEFILE ON)
SET(CMAKE_C_COMPILER mpicc)
ADD_DEFINITIONS(-g)

INCLUDE(CMakeFortran.txt)

ADD_DEFINITIONS(-DMANGLE_F77_${SCHEME})

ADD_LIBRARY (nbc ${LibNBC_SRC} ${LibNBC_COLLALGS_SRC})

