#!/bin/bash
#gfortran -Wall -Wextra -Wtabs -fPIC -g -fcheck=all -fbacktrace -ffree-line-length-0 cylinder_system.f08 dranxor2.f95 -o cylinder_system.exe #use this one for development
gfortran -Wall -Wextra -Wtabs -fPIC -O3 -march=native -ffast-math -funroll-loops -ffree-line-length-0 cylinder_system.f08 dranxor2.f95 -o cylinder_system_final.exe #use this one for release
