program octopus_mpi
  use loct_oct_m

  write(*, *) "Error: The 'octopus_mpi' command is now obsolete,"
  write(*, *) "       please use the 'octopus' command."

  call loct_exit_failure()

end program octopus_mpi
