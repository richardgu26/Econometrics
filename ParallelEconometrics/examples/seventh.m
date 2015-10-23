# shows use of MPI_Iprobe
LAM_Init(1);
global TAG NEWORLD;
# [info, size] = MPI_Comm_size(NEWORLD);
# size


# send objects and the command to run
# note how the node will pause a while before sending the result
cmd = 'pause(10); a = rand(3, 3); b = inv(a); MPI_Send(b, 0, TAG, NEWORLD);';
NumCmds_Send({"cmd"}, {cmd});

b = eye(3); # need a proper sized container to hold what's received

done = false;
while (!done)
	[info ready] = MPI_Iprobe(1, TAG, NEWORLD); # check if message pending
	if ready
		printf("OK, now I'm ready\n");
		MPI_Recv(b, 1, TAG, NEWORLD);
		done = true;
	else
		printf("Hold on, I'm not ready yet!\n");
		pause(1);
	endif
endwhile
printf("b received from rank 1:\n");
b
LAM_Finalize;
