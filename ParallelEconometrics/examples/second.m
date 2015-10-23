MPI_Init;
# use help to see what MPI_COMM_WORLD is
[info, size] = MPI_Comm_size(MPI_COMM_WORLD);
size
[info, rank] = MPI_Comm_rank(MPI_COMM_WORLD);
rank
MPI_Finalize;
MPITB_Finalize; # this is needed to avoid confusing Octave
