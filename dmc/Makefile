FC = gfortran
FFLAGS =
SRC = dmc.f95 normrand.f95 system.f95
OBJ = ${SRC:.f95=.o}

%.o: %.f95
	$(FC) $(FFLAGS) -o $@ -c $<

dmc: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)
  
clean:
	-rm -f *.o
