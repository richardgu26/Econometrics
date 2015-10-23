# on each rank, this adds the rank to a constant and sends back to rank 0
# play around increasing size of communicator
nodes = 3;
LAM_Init(nodes);
global TAG, NEWORLD;
# [info, size] = MPI_Comm_size(NEWORLD);
# size


dim = 3;
# send objects and the command to run on each rank
cmd = 'a = rand(dim, dim); b = inv(a); MPI_Send(b, 0, TAG, NEWORLD);';
NumCmds_Send({"dim", "cmd"}, {dim, cmd});

for i = 1:nodes
b = eye(dim); # need a proper sized container to hold what's received
MPI_Recv(b, i, TAG, NEWORLD);
printf("b received from rank %d:\n", i);
b
endfor

LAM_Finalize;
