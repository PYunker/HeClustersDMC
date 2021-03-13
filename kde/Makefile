SRCS=$(wildcard *.f90)
OBJS=$(patsubst %.f90,%.o,$(SRCS))

FC = mpifort
FLFLAGS = -llapack -lblas
FCFLAGS = -c -Wall -Wextra -Wconversion -fcheck=bounds -ffpe-trap=invalid -fopenmp
PROGRAM = mse

default : $(PROGRAM)

$(OBJS) : %.o : %.f90
    $(FC) $(FCFLAGS) -o $@ $<

$(PROGRAM) : $(OBJS)
    $(FC) $(FLFLAGS) -o $@ $^

debug:
    @echo "SRCS = $(SRCS)"
    @echo "OBJS = $(OBJS)"
    @echo "PROGRAM = $(PROGRAM)"

clean:
    rm -rf $(OBJS) $(PROGRAM) *.mod

.PHONY: debug default clean

# Main program depends on all modules
$(PRG_OBJ) : $(MOD_OBJS)


mod_blocks.o mod_allocations.o : mod_shared.o
