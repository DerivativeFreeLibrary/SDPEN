-----------------------------------------------------------
 How to use the derivative-free optimizer SDPEN (MATLAB version)
-----------------------------------------------------------

1- Gunzip and untar the archive seqpen_matlab.tar.gz in a folder on your computer.

The folder contains the following files:

- sdpen.m		(optimizer routine)
- example.m		(main script to run an example problem)
- fob_con.m 	(function used in the example, it computes f.obj
                 and constraint function values)

-----------------------------------------------------------
 How to use the derivative-free optimizer SDPEN (F90 version)
-----------------------------------------------------------

1- Gunzip and untar the archive seqpen_f90.tar.gz in a folder on your computer.

2- Edit file problem.f90 to define your own problem.
   In particular, you need to define the routines
   setdim 	which sets the dimensions of the problem
   setbounds 	which sets upper and lower bounds on the variables
   startp	 which sets the starting point
   funob 	which returns the obj. function value
   fconstr 	which returns the constraint values

3- At command prompt execute 

     $> make
 
   which will create the executable 'penseq'

4- execute

     $> ./penseq

-----------------------------------------------------------
 How to use the derivative-free optimizer SDPEN (C version)
-----------------------------------------------------------

1- Gunzip and untar the archive seqpen_c.tar.gz in a folder on your computer.

2- Edit file problem.c to define your own problem.
   In particular, you need to change the routines
   setdim 	which sets the dimensions of the problem
   setbounds 	which sets upper and lower bounds on the variables
   startp	which sets the starting point
   fobj 	which returns the obj. function value
   fconstr 	which returns the constraint values

3- At command prompt execute 

     $> make
 
   which will create the executable 'seqpen'

4- execute

     $> ./seqpen

Notes on the C version.

Two ways of allocation are included for the vectors.
A static allocation, given in seqpen.h, and a dynamic allocation also defined
in seqpen.h and carried out in seqpenmain.c. Appropriate sections must
be commented out or activated, as described in seqpen.h and
seqpenmain.c. Currently, static allocation has been specified,
with max sizes specified in basic.h.

