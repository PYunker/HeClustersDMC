# This is admittedly quite messy ($(OBJDIR) spammed everywhere etc...)

SRCDIR = src
OBJDIR = obj
BINDIR = bin

SRCS := $(wildcard $(SRCDIR)/*.f90)
# OBJS := $(patsubst %.f90,%.o,$(SRCS))
OBJS := $(SRCS:$(SRCDIR)/%.f90=$(OBJDIR)/%.o)

FC = mpifort
FLFLAGS = -fopenmp
FCFLAGS = -fopenmp -fintrinsic-modules-path $(OBJDIR) -J$(OBJDIR)
PROGRAM = mse

default : $(BINDIR)/$(PROGRAM)

$(OBJS) : $(OBJDIR)/%.o : $(SRCDIR)/%.f90
	$(FC) $(FCFLAGS) -o $@ -c $<

$(OBJDIR)/main.o : $(OBJDIR)/kde.o $(OBJDIR)/distributions.o $(OBJDIR)/grid_utils.o $(OBJDIR)/metropolis.o
$(OBJDIR)/metropolis.o : $(OBJDIR)/distributions.o

$(BINDIR)/$(PROGRAM) : $(OBJS)
	$(FC) $(FLFLAGS) -o $@ $^ -llapack -lblas
	@echo "linking complete"

debug:
	@echo "SRCS = $(SRCS)"
	@echo "OBJS = $(OBJS)"
	@echo "PROGRAM = $(PROGRAM)"

clean:
	rm -rf $(OBJS) $(PROGRAM) $(OBJDIR)/*.mod

remove: clean
	rm -f $(BINDIR)/$(PROGRAM)

.PHONY: debug default clean remove
