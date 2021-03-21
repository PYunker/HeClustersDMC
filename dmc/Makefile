FC = gfortran-10
FFLAGS = -fsanitize=leak
SRC = dmc_utils.f90  main.f90
OBJ = ${SRC:.f90=.o}

%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $<

dmc: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ) -llapack -lblas

main.o: dmc_utils.o

debug:
	@echo $(OBJ)
	@echo $(SRC)
  
clean:
	-rm -f *.o *.mod dmc populace_dmc
