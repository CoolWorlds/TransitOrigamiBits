#!/bin/bash
#

rm -f f*/plotit_chi_txt.txt
rootfile='../plotit_chi_txt.txt'

head -3 $rootfile > f0.1/plotit_chi_txt.txt
head -3 $rootfile > f0.2/plotit_chi_txt.txt
head -3 $rootfile > f0.3/plotit_chi_txt.txt
head -3 $rootfile > f0.4/plotit_chi_txt.txt
head -3 $rootfile > f0.5/plotit_chi_txt.txt
head -3 $rootfile > f0.6/plotit_chi_txt.txt
head -3 $rootfile > f0.7/plotit_chi_txt.txt
head -3 $rootfile > f0.8/plotit_chi_txt.txt
head -3 $rootfile > f0.9/plotit_chi_txt.txt
head -3 $rootfile > f1.0/plotit_chi_txt.txt

thisfolder=`pwd`

line=`echo "thisfolder=\""$thisfolder"/f0.1/\""`
echo $line >> f0.1/plotit_chi_txt.txt
line=`echo "thisfolder=\""$thisfolder"/f0.2/\""`
echo $line >> f0.2/plotit_chi_txt.txt
line=`echo "thisfolder=\""$thisfolder"/f0.3/\""`
echo $line >> f0.3/plotit_chi_txt.txt
line=`echo "thisfolder=\""$thisfolder"/f0.4/\""`
echo $line >> f0.4/plotit_chi_txt.txt
line=`echo "thisfolder=\""$thisfolder"/f0.5/\""`
echo $line >> f0.5/plotit_chi_txt.txt
line=`echo "thisfolder=\""$thisfolder"/f0.6/\""`
echo $line >> f0.6/plotit_chi_txt.txt
line=`echo "thisfolder=\""$thisfolder"/f0.7/\""`
echo $line >> f0.7/plotit_chi_txt.txt
line=`echo "thisfolder=\""$thisfolder"/f0.8/\""`
echo $line >> f0.8/plotit_chi_txt.txt
line=`echo "thisfolder=\""$thisfolder"/f0.9/\""`
echo $line >> f0.9/plotit_chi_txt.txt
line=`echo "thisfolder=\""$thisfolder"/f1.0/\""`
echo $line >> f1.0/plotit_chi_txt.txt

tail -70 $rootfile >> f0.1/plotit_chi_txt.txt
tail -70 $rootfile >> f0.2/plotit_chi_txt.txt
tail -70 $rootfile >> f0.3/plotit_chi_txt.txt
tail -70 $rootfile >> f0.4/plotit_chi_txt.txt
tail -70 $rootfile >> f0.5/plotit_chi_txt.txt
tail -70 $rootfile >> f0.6/plotit_chi_txt.txt
tail -70 $rootfile >> f0.7/plotit_chi_txt.txt
tail -70 $rootfile >> f0.8/plotit_chi_txt.txt
tail -70 $rootfile >> f0.9/plotit_chi_txt.txt
tail -70 $rootfile >> f1.0/plotit_chi_txt.txt