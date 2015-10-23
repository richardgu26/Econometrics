# remember to use "help" to learn about the functions
info = MPI_Init;
info
[info, flag] = MPI_Initialized;
flag
MPI_Finalize;
MPITB_Finalize; # this is needed to avoid confusing Octave
[info, flag] = MPI_Initialized;
flag
