#!/bin/bash
gfortran -Wall -Wextra -Wtabs -mcmodel=large -fPIC -g -fcheck=all -fbacktrace -ffree-line-length-none cylinder_system.f08 dranxor2.f95 -o cylinder_system_debug.exe #use this one for development
echo ""
echo "##########################"
echo "Debugging binary compiled."
echo "##########################"
echo ""
gfortran -Wall -Wextra -Wtabs -mcmodel=large -fPIC -O3 -march=native -ffast-math -funroll-loops -ffree-line-length-none cylinder_system.f08 dranxor2.f95 -o cylinder_system.exe #use this one for release
echo ""
echo "########################"
echo "Release binary compiled."
echo "########################"
echo ""
#gfortran -Wall -Wextra -Wtabs -mcmodel=large -fPIC -g -fcheck=all -fbacktrace -fopenmp -ffree-line-length-none cylinder_system.f08 dranxor2.f95 -o cylinder_system_parallel_debug.exe #use this one for parallel development
#echo ""
#echo "####################################"
#echo "Debugging parallel binary compiled."
#echo "####################################"
#echo ""
#gfortran -Wall -Wextra -Wtabs -mcmodel=large -fPIC -O3 -march=native -ffast-math -funroll-loops -fopenmp -ffree-line-length-none cylinder_system.f08 dranxor2.f95 -o cylinder_parallel_system.exe #use this one for parallel release
#echo ""
#echo "#################################"
#echo "Release parallel binary compiled."
#echo "#################################"
#echo ""
