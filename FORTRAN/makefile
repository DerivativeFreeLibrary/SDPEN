FC = gfortran
RM = rm -f

FFLAGS = -O3

OBJS = problem.o SDPEN.o PenSeq-main.o

all:  $(OBJS) 
	$(FC) -o penseq $(OBJS)

.SUFFIXES : .f90 .o

.f90.o: $* ; $(FC) $(FFLAGS) -c $*.f90

clean: 
	$(RM) *.mod
	$(RM) *.o
	$(RM) penseq

