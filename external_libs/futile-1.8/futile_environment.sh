prefix=/local/binaries/gfortran-dist-test/build/install
exec_prefix=${prefix}
. "${exec_prefix}/bin/unienv.sh"
unienv LD_LIBRARY_PATH ${exec_prefix}/lib
unienv PYTHONPATH ${exec_prefix}/lib/python2.7/site-packages
