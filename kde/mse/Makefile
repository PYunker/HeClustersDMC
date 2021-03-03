FC = mpif90
FFLAGS =# -ffpe-trap=invalid -Wall
SRC = metropolis.f90 grid_utils.f90 kde.f90 main.f90
OBJ = ${SRC:.f90=.o}

%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $<

mse: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ) -llapack -lblas
  
.PHONY clean:
	-rm -f *.o *.mod mse 
