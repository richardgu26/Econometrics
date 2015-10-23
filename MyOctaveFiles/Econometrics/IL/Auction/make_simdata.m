# this is to solve a model many times using different parameter values
# drawn from a prior distribution. This runs on a cluster. The frontend
# node (MPI rank 0) does not participate in the solution, it only collects
# the results and writes to disk. You should run this using something like
# mpirun -np X --hostfile /foo/bar octave -q --eval make_simdata
# where X is the total number of cores in the cluster, PLUS ONE. That way,
# all cores will be used in computations, and one will have the added duty
# of serving as the frontend

verbose = 1;

% initial parameter bounds: max mean value is 150, whatever the quality
lb_param_ub = [
-1     0.5    3	% beta0
0.0    0.5   2	% beta1
];

if not(MPI_Initialized) MPI_Init; endif
CW = MPI_Comm_Load("NEWORLD");
node = MPI_Comm_rank(CW);
nodes = MPI_Comm_size(CW);
mytag = 48;

outsamp_per_node = ceil(outsamp/nodes);


model_params0 = lb_param_ub(:,2);
lb = lb_param_ub(:,1);
ub = lb_param_ub(:,3);


###############################################################################
############## you do not need to alter anything below this line ##############
###############################################################################


if node
	more_please = 1;
	while more_please
		# break it up to do intermediate writes
		for j = 1:n_pooled
			if outsamp
				model_params = model_params0;
			else	
				model_params = rand(rows(ub),1).*(ub-lb) + lb;
			endif

			% the model
			theta1 = model_params(1,:);
			theta2 = model_params(2,:);
			N = 6;
			# quality of good
			x = rand(n,1);
			# valuations drawn from exponetial mean phi
			phi = exp(theta1 + theta2*x);
			# highest valuation
			r = rand(n,6);
			r = min(r')';
			v = -log(r).*phi;
			# get winning bid
			z = v./phi;
			D = exp(-5*z).*(60*exp(5*z) + 300*phi .* exp(4*z) - 300*phi .* exp(3*z) ...
			+ 200*phi .* exp(2*z) - 75*phi .* exp(z) + 12*phi)/60 - 137*phi/60;
			b = v - D ./ ((1 - exp(-v./phi)).^(N-1));
			b = b.*(b>0);
			% the aux stat
			b = 0.01*(b<0.01) + b.*(b>0.01);
			y = log(b);
			z = [ones(n,1) x];
			bhat = z\y;
			e = y - z*bhat;
			s = log(e'*e/(n-2));
			m1 = mean(log(b));
			m2 = std(log(b));
			m3 = mean((log(b)-m1).^3);
			contrib = [outsamp theta1 theta2 bhat' s m1 m2 m3];
			if (i==1) && (j==1) contribs = zeros(n_pooled, columns(contrib)); endif
			contribs(j,:) = contrib;
		endfor
		# check if we're done
		if (MPI_Iprobe(0, mytag+1, CW))
			junk = MPI_Recv(0, mytag+1, CW);
			more_please = false;
		else
			MPI_Send(contribs, 0, mytag, CW);
		endif		
	endwhile
else # frontend
	received = 0;
	rescaled = false;
	done = false;
	while received < reps
		# retrieve results from compute nodes
		pause(0.05);
		for i = 1:nodes-1
			# compute nodes have results yet?
			ready = false;
			ready = MPI_Iprobe(i, mytag, CW); # check if message pending
			if ready
				# get it if it's there
				contribs = MPI_Recv(i, mytag, CW);
			       	need = reps - received;
				received = received + n_pooled;
				# truncate?
				if n_pooled  >= need
						contribs = contribs(1:need,:);
						done = true;
				endif
				# write to output file
				FN = fopen (outfile, "a");
				if (FN < 0) error ("make_simdata: couldn't open output file %s", outfile); endif
				for j = 1:rows(contribs)
					fprintf(FN, "%f ", contribs(j,:));
					fprintf(FN, "\n");
				endfor
				fclose(FN);
				system('sync');
				if verbose
					printf("\nContribution received from node%d.  Received so far: %d\n", i, received);
				endif
			endif
		endfor
		if done
			# tell compute nodes to stop loop
			for i = 1:(nodes-1)
				MPI_Send("stop!",i, mytag+1,CW); # send out message to stop
				ready = MPI_Iprobe(i, mytag, CW); # get last messages
				if ready contribs = MPI_Recv(i, mytag, CW); endif
			endfor
			break;
		endif
	endwhile
endif
if not(MPI_Finalized) MPI_Finalize; endif

