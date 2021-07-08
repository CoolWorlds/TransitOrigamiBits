#!/bin/bash
#

math='/Applications/Mathematica.app/Contents/MacOS/wolframscript'

cd f0.1/
$math -script plotit_chi_txt.txt
echo 'f0.1 complete'
cd ../f0.2/
$math -script plotit_chi_txt.txt
echo 'f0.2 complete'
cd ../f0.3/
$math -script plotit_chi_txt.txt
echo 'f0.3 complete'
cd ../f0.4/
$math -script plotit_chi_txt.txt
echo 'f0.4 complete'
cd ../f0.5/
$math -script plotit_chi_txt.txt
echo 'f0.5 complete'
cd ../f0.6/
$math -script plotit_chi_txt.txt
echo 'f0.6 complete'
cd ../f0.7/
$math -script plotit_chi_txt.txt
echo 'f0.7 complete'
cd ../f0.8/
$math -script plotit_chi_txt.txt
echo 'f0.8 complete'
cd ../f0.9/
$math -script plotit_chi_txt.txt
echo 'f0.9 complete'
cd ../f1.0/
$math -script plotit_chi_txt.txt
echo 'f1.0 complete'
cd ..

echo 'DONE2'