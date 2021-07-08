#!/bin/bash
#

fortran='gfortran'

# move over program files
rm -f *.o *.mod *.pdf PRI_full.jam PRI_fold.jam *.pdf

$fortran -O3 -c params.f90
$fortran -O3 -c jasmine.f90
$fortran -O3 -c mandel.f
$fortran -O3 -c luna.f90
$fortran -O3 -c transit.f90
$fortran -O3 -o plotit plotit.f90 *.o

./plotit
mv PRI_full.jam main_output.jam

$fortran -O3 -c params.f90
$fortran -O3 -c jasmine.f90
$fortran -O3 -c mandel.f
$fortran -O3 -c luna.f90
$fortran -O3 -c transit.f90
$fortran -O3 -o plotit2 plotit2.f90 *.o

./plotit2
mv PRI_full.jam ghost_output.jam