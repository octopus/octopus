! This file has been processed by configure, do not change it.
! The source file is ${BigDFT sources}/src/modules/configure.inc.in.

!Memory limit value in GB. It stops EVERYTHING if some process passes such limit
!For no memory limit, leave it to zero
real(kind=4), parameter :: memorylimit = 0

!Convenience variables.
character(len = *), parameter :: package_name      = "BigDFT - DFT over wavelets"
character(len = *), parameter :: package_tarname   = "bigdft"
character(len = *), parameter :: package_version   = "1.7.6"
character(len = *), parameter :: package_string    = "BigDFT - DFT over wavelets 1.7.6"
character(len = *), parameter :: package_bugreport = "Damien.Caliste@cea.fr"

