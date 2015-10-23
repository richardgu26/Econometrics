# play around increasing size of communicator
LAM_Init(15);
global NEWORLD;
[info, size] = MPI_Comm_size(NEWORLD);
size

# rank of compute node
cmd = '[info, rank] = MPI_Comm_rank(NEWORLD); printf("my rank is %d\n", rank)';
NumCmds_Send({"cmd"},{cmd});

# rank of frontend node
eval(cmd);
LAM_Finalize;
