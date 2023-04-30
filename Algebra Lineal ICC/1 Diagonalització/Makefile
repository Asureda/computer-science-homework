##Makefile per una simulacio de Dinamica Molecular


# Compilador
F90=mpifort

# Main program
TARGET=main

# Optimitzador
OPT=-O

# Flags
FLAGS=-Wall  -fbounds-check

# Grafics d'energia, pressio, temperatura i funcio de distribució radial
#energy.png : termodynamics_reduced.dat termodynamics_real.dat distrib_func.dat
#    @echo "Generants plots amb els resultats..."
#    gnuplot red.gnu
#    gnuplot real.gnu
#    gnuplot gr.gnu

# Execucio del main program
#termodynamics_reduced.dat : $(TARGET) config.dat constants.dat parameters.dat
#    @echo "Executant el programa amb els valors de l'input..."
#    ./$(TARGET)

# Compilacio del main program
$(TARGET) : diag.o $(TARGET).o
	@echo "From MAKEFILE. Compiling main program"
	$(F90) $(OPT) $(FLAGS) -o $(TARGET) $(TARGET).o diag.o

# Tots els fitxers amb extensio .f90 son compilats a objectes .o
%.o : %.f90
	@echo "From MAKEFILE. Compiling all modules"
	$(F90) $(OPT) $(FLAGS) -c $<


## help: instruccions sobre l'us d'aquest Makefile
.PHONY : help
help :
	@sed -n 's/^##//p' Makefile


## backup : fer una copia comprimida del codi basic
.PHONY : backup
backup:
	tar -czvf "backup.tar.gz" *.f90

## clean_all : regla per fer un clean dels objectes executables, resultats i imatges
.PHONY : clean_all
clean_all:
	@echo Removing compiled objects and results
	rm *.o *.dat *.eps *.xyz *.x*.mod

## clean_exe : regla per fer un clean nomes dels objectes executables
.PHONY : clean_exe
clean_exe:
	@echo Removing compiled objects only
	rm *.o *.mod
