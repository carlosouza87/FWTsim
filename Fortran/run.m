clear all;close all;clc
dos('gfortran sysvar.f90 lib_array.f90 syssub.f90 FWTsim.f90 -o FWTsim -Og');
dos('FWTsim');
plt