#!/bin/bash
gfortran -Wextra -fPIC -g -fcheck=all -fbacktrace -ffree-line-length-0 cylinder_system.f08 dranxor2.f95 -o cylinder_system.exe #development
#gfortran -Wextra -fPIC -Werror -O3 -march=native -ffast-math -funroll-loops -ffree-line-length-0 cylinder_system.f08 dranxor2.f95 -o cylinder_system_final.exe #final
