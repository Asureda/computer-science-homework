opt: linear_interpolation.f90
	gfortran -o Ex3_4_1 Ex3_4_1.f90 -O3

check: linear_interpolation.f90
	gfortran -Wall -Wextra -fbounds-check -o Ex3_4_1 Ex3_4_1.f90

backup: Ex3_4_1.f90
	rm -rf ./backup
	mkdir ./backup
	cp Ex3_4_1 ./backup/
	cp Ex3_4_1 ./backup/


