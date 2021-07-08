#!/bin/bash
#

fortran='gfortran'

# move over program files
rm -f */*.pdf */*.o */*.mod */plotit */plotit2
rm -f */PRI_fold.jam */main_output.jam */ghost_output.jam
rm -f */chi2_scales.jam */Msp_grid.dat */ttv.dat

cd f0.1/
./runme.sh
cd ../f0.2/
./runme.sh
cd ../f0.3/
./runme.sh
cd ../f0.4/
./runme.sh
cd ../f0.5/
./runme.sh
cd ../f0.6/
./runme.sh
cd ../f0.7/
./runme.sh
cd ../f0.8/
./runme.sh
cd ../f0.9/
./runme.sh
cd ../f1.0/
./runme.sh
cd ..

echo 'DONE'