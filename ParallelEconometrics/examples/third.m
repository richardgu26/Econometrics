LAM_Init(1);
global NEWORLD;
[info, size] = MPI_Comm_size(NEWORLD);
size
[info, rank] = MPI_Comm_rank(NEWORLD);
rank
LAM_Finalize;
