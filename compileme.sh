#!/bin/bash

#gfortran -fopenmp mod_global_variables.f90 mod_read.f90 mod_utils.f90 mod_rates.f90 MAIN.f90 -o analysis.exe
pgf95 mod_global_variables.f90 mod_read.f90 mod_utils.f90 mod_rates.f90 mod_density.f90 mod_misc.f90 MAIN.f90 \
-Minfo -Mneginfo -time -fast -O3 -Miomutex -Mconcur=allcores -mp=allcores -Munroll -Mvect -llapack -lblas -o analysis.exe
