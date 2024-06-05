gfortran -static -fopenmp -o3 -o PeriCrush_v10.exe ../src/PeriCrush_v10.f90
del *.mod
PeriCrush_v10.exe