if not(MPI_Initialized) MPI_Init; endif

S = 500; # replications for II
n = 80; # sample size
outfile = "AuctionII.80new";
wrapper_args = {n, S};
reps = 5e3;
reps = 5000;
n_pooled = 1;
montecarlo("wrapper", wrapper_args, reps, outfile, n_pooled, true, true);
if not(MPI_Finalized) MPI_Finalize; endif

