FC = gfortran
FFLAGS =
SRC = ../dmc.f95 ../normrand.f95 har_os_nd.f95
OBJ = ${SRC:.f95=.o}

%.o: %.f95
	$(FC) $(FFLAGS) -o $@ -c $<

dmc: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)
  
clean:
	-rm -f ../*.o *.o *.mod dmc populace_dmc p_tran
