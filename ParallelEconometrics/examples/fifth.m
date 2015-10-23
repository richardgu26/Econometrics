# on each rank, this adds the rank to a constant and sends back to rank 0
# play around increasing size of communicator
nodes = 3;
LAM_Init(nodes);
global TAG NEWORLD;
[info, size] = MPI_Comm_size(NEWORLD);
size

a = 1;

# send a and a command to run on each rank
cmd = '[info, rank] = MPI_Comm_rank(NEWORLD); a = a + rank; MPI_Send(a, 0, TAG, NEWORLD);';
NumCmds_Send({"a", "cmd"},{a, cmd});

for i = 1:nodes
MPI_Recv(a, i, TAG, NEWORLD);
printf("value of a received from rank %d: %d\n", i, a);
endfor

LAM_Finalize;
