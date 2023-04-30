#!/bin/bash

#Compliler for C and Fortran 
comp1='gcc'
comp2='gfortran'

# Optimization Flags



# Warning and Error flags gfortran
gflags_hard='-Wno-tabs -Wall -Wextra -Warray-temporaries -fbounds-check -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan'
gflags_soft='-Wall -fbounds-check -Wno-tabs'
# Warning and Error flags intel compiler
iflags_hard='-check all -fpe0 -warn -traceback -debug extended'
iflags_soft=''

flags=
#$gflags_hard
#opt=
opt='-O'

$comp2 $opt -c $flags schrodinger.f
$comp2 schrodinger.o $opt $flags -o r_main
rm *.o
rm *.mod