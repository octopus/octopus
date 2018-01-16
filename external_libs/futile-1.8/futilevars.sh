echo "While sourcing this script the following environment variable will be set:"
prefix=/local/binaries/gfortran-dist-test/build/install
exec_prefix=${prefix}
sourcefiles="${exec_prefix}/bin/futile_environment.sh"
for file in $sourcefiles
do
sh $file
eval `sh $file`
done
$(return >/dev/null 2>&1)
if [ "$?" -eq "0" ]
then
    echo "...done."
else
    echo "ERROR: This script is NOT sourced! Execute 'source" $0 "' instead"
fi
