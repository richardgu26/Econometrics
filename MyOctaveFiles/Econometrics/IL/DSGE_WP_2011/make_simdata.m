# this is to solve a model many times using different parameter values
# drawn from a prior distribution. This runs on a cluster. The frontend
# node (MPI rank 0) does not participate in the solution, it only collects
# the results and writes to disk. You should run this using something like
# mpirun -np X --hostfile /foo/bar octave -q --eval make_simdata
# where X is the total number of cores in the cluster, PLUS ONE. That way,
# all cores will be used in computations, and one will have the added duty
# of serving as the frontend

# the path of output file needs to be specified, because I'm running from a ramdisk.
# I want the output to go to permanent storage
outfile = "simdata.uniform";
reps = 100;			# total number of simulations
percentage_at_design = 0.2;  	# % of total simulations at design point, rest from param space
n_pooled = 1;
verbose = 1;
if not(MPI_Initialized) MPI_Init; endif
CW = MPI_Comm_Load("NEWORLD");
node = MPI_Comm_rank(CW);
nodes = MPI_Comm_size(CW);
mytag = 48;

%----------------------------------------------------------------
% 2. Calibration
%----------------------------------------------------------------
alpha   = 0.33;
beta    = 0.99;
delta   = 0.023;
psi     = 1.75;
rho     = 0.95;
sigma   = 0.010448;
epsilon = 10.0;
model_params1 = [alpha; beta; delta; psi; rho; sigma; epsilon];

first_time = true;
if node
	more_please = 1;
	while more_please
		# break it up to do intermediate writes
		for j = 1:n_pooled
			# re-establish the state of the RNG, dynare seems to synchronize this across MPI nodes
			if !first_time
				if     node==1
					load state1;
					rand("state", s);
				elseif node==2
					load state2;
					rand("state", s);
				elseif node==3
					load state3;
					rand("state", s);
				elseif node==4
					load state4;
					rand("state", s);
				elseif node==5
					load state5;
					rand("state", s);
				elseif node==6
					load state6;
					rand("state", s);
				elseif node==7
					load state7;
					rand("state", s);
				elseif node==8
					load state8;
					rand("state", s);
				elseif node==9
					load state9;
					rand("state", s);
				elseif node==10
					load state10;
					rand("state", s);
				elseif node==11
					load state11;
					rand("state", s);
				elseif node==12
					load state12;
					rand("state", s);
				elseif node==13
					load state13;
					rand("state", s);
				elseif node==14
					load state14;
					rand("state", s);
				elseif node==15
					load state15;
					rand("state", s);
				elseif node==16
					load state16;
					rand("state", s);
				elseif node==17
					load state17;
					rand("state", s);
				elseif node==18
					load state18;
					rand("state", s);
				elseif node==19
					load state19;
					rand("state", s);
				elseif node==20
					load state20;
					rand("state", s);
				elseif node==21
					load state21;
					rand("state", s);
				elseif node==22
					load state22;
					rand("state", s);
				elseif node==23
					load state23;
					rand("state", s);
				elseif node==24
					load state24;
					rand("state", s);
				elseif node==25
					load state25;
					rand("state", s);
				elseif node==26
					load state26;
					rand("state", s);
				elseif node==27
					load state27;
					rand("state", s);
				elseif node==28
					load state28;
					rand("state", s);
				elseif node==29
					load state29;
					rand("state", s);
				elseif node==30
					load state30;
					rand("state", s);
				elseif node==31
					load state31;
					rand("state", s);
				elseif node==32
					load state32;
					rand("state", s);
				endif
			endif

			first_time=false; # from now on we want to keep the RNG sequence moving along properly

			# most draws will be from param space, but some from true values
			a = rand(1,1);
			if (a < percentage_at_design)
				model_params = model_params1;
				insamp = 1;
			else
				# these are fairly informative priors
				# the means and st. dev. that result are noted
				#alpha   = betarnd(1.9871e+02 , 3.6904e+02     ); # m=0.35 s=0.02
				#beta    = betarnd(1.9110e+02 , 3.9000e+00     ); # m=0.98 s=0.01
				#delta   = betarnd(2.4350e+01 , 9.4965e+02     ); # m=0.025 s=0.005
				#psi     = gamrnd (3.0625e+02 , 5.7143e-03     ); # m=1.75 s=0.1
				#rho     = betarnd(1.1186e+02 , 5.8875e+00     ); # m=0.95 s=0.02
				#epsilon = gamrnd (1.0000e+04 , 1.0000e-03     ); # m=10 s=0.1
				#siginv =  gamrnd (3.0000e+00 , 5.0000e+01     );
				#sigma = 1/siginv;				 # m = 0.01 s=0.01
				#sigma = 0.015*rand(1,1) + 0.005; # uniform 0.005-0.02
				#model_params = [alpha; beta; delta; psi; rho; sigma; epsilon];
				# bounds of parameter space for uninformative (uniform) prior
				lb = [0.15; 0.95; 0.005; 1; 0.85; 0.005; 9];
				ub = [0.4; 0.999; 0.06; 3; 0.99; 0.04; 13];
				model_params = rand(7,1).*(ub-lb) + lb;
				insamp = 0;
			endif

			# break into pieces and save
			alpha = model_params(1,:);
			beta = model_params(2,:);
			delta = model_params(3,:);
			psi = model_params(4,:);
			rho = model_params(5,:);
			sigma = model_params(6,:);
			epsilon = model_params(7,:);


			# there are a number of .mod files that each read a different parameter file
			# this writes out the different parameter files so that a number of solutions
			# can be done in parallel
			s = rand("state");
			if     node==1
				save parameterfile1  alpha beta delta psi rho sigma epsilon;
				save state1 s;
			elseif node==2
				save parameterfile2  alpha beta delta psi rho sigma epsilon;
				save state2 s;
			elseif node==3
				save parameterfile3  alpha beta delta psi rho sigma epsilon;
				save state3 s;
			elseif node==4
				save parameterfile4  alpha beta delta psi rho sigma epsilon;
				save state4 s;
			elseif node==5
				save parameterfile5  alpha beta delta psi rho sigma epsilon;
				save state5 s;
			elseif node==6
				save parameterfile6  alpha beta delta psi rho sigma epsilon;
				save state6 s;
			elseif node==7
				save parameterfile7  alpha beta delta psi rho sigma epsilon;
				save state7 s;
			elseif node==8
				save parameterfile8  alpha beta delta psi rho sigma epsilon;
				save state8 s;
			elseif node==9
				save parameterfile9  alpha beta delta psi rho sigma epsilon;
				save state9 s;
			elseif node==10
				save parameterfile10 alpha beta delta psi rho sigma epsilon;
				save state10 s;
			elseif node==11
				save parameterfile11 alpha beta delta psi rho sigma epsilon;
				save state11 s;
			elseif node==12
				save parameterfile12 alpha beta delta psi rho sigma epsilon;
				save state12 s;
			elseif node==13
				save parameterfile13 alpha beta delta psi rho sigma epsilon;
				save state13 s;
			elseif node==14
				save parameterfile14 alpha beta delta psi rho sigma epsilon;
				save state14 s;
			elseif node==15
				save parameterfile15 alpha beta delta psi rho sigma epsilon;
				save state15 s;
			elseif node==16
				save parameterfile16 alpha beta delta psi rho sigma epsilon;
				save state16 s;
			elseif node==17
				save parameterfile17 alpha beta delta psi rho sigma epsilon;
				save state17 s;
			elseif node==18
				save parameterfile18 alpha beta delta psi rho sigma epsilon;
				save state18 s;
			elseif node==19
				save parameterfile19 alpha beta delta psi rho sigma epsilon;
				save state19 s;
			elseif node==20
				save parameterfile20 alpha beta delta psi rho sigma epsilon;
				save state20 s;
			elseif node==21
				save parameterfile21 alpha beta delta psi rho sigma epsilon;
				save state21 s;
			elseif node==22
				save parameterfile22 alpha beta delta psi rho sigma epsilon;
				save state22 s;
			elseif node==23
				save parameterfile23 alpha beta delta psi rho sigma epsilon;
				save state23 s;
			elseif node==24
				save parameterfile24 alpha beta delta psi rho sigma epsilon;
				save state24 s;
			elseif node==25
				save parameterfile25 alpha beta delta psi rho sigma epsilon;
				save state25 s;
			elseif node==26
				save parameterfile26 alpha beta delta psi rho sigma epsilon;
				save state26 s;
			elseif node==27
				save parameterfile27 alpha beta delta psi rho sigma epsilon;
				save state27 s;
			elseif node==28
				save parameterfile28 alpha beta delta psi rho sigma epsilon;
				save state28 s;
			elseif node==29
				save parameterfile29 alpha beta delta psi rho sigma epsilon;
				save state29 s;
			elseif node==30
				save parameterfile30 alpha beta delta psi rho sigma epsilon;
				save state30 s;
			elseif node==31
				save parameterfile31 alpha beta delta psi rho sigma epsilon;
				save state31 s;
			elseif node==32
				save parameterfile32 alpha beta delta psi rho sigma epsilon;
				save state32 s;
			endif

			# solve model, using the appropriate parameter file
			if node==1
				dynare RBC1 noclearall;
			elseif node==2
				dynare RBC2 noclearall;
			elseif node==3
				dynare RBC3 noclearall;
			elseif node==4
				dynare RBC4 noclearall;
			elseif node==5
				dynare RBC5 noclearall;
			elseif node==6
				dynare RBC6 noclearall;
			elseif node==7
				dynare RBC7 noclearall;
			elseif node==8
				dynare RBC8 noclearall;
			elseif node==9
				dynare RBC9 noclearall;
			elseif node==10
				dynare RBC10 noclearall;
			elseif node==11
				dynare RBC11 noclearall;
			elseif node==12
				dynare RBC12 noclearall;
			elseif node==13
				dynare RBC13 noclearall;
			elseif node==14
				dynare RBC14 noclearall;
			elseif node==15
				dynare RBC15 noclearall;
			elseif node==16
				dynare RBC16 noclearall;
			elseif node==17
				dynare RBC17 noclearall;
			elseif node==18
				dynare RBC18 noclearall;
			elseif node==19
				dynare RBC19 noclearall;
			elseif node==20
				dynare RBC20 noclearall;
			elseif node==21
				dynare RBC21 noclearall;
			elseif node==22
				dynare RBC22 noclearall;
			elseif node==23
				dynare RBC23 noclearall;
			elseif node==24
				dynare RBC24 noclearall;
			elseif node==25
				dynare RBC25 noclearall;
			elseif node==26
				dynare RBC26 noclearall;
			elseif node==27
				dynare RBC27 noclearall;
			elseif node==28
				dynare RBC28 noclearall;
			elseif node==29
				dynare RBC29 noclearall;
			elseif node==30
				dynare RBC30 noclearall;
			elseif node==31
				dynare RBC31 noclearall;
			elseif node==32
				dynare RBC32 noclearall;
			endif

			# get a simulation of length 80 (20 years quarterly), and compute aux. statistic
			data = [y l c k w r] ;
			data = data(101:180,:);
			Z = [aux_stat(data)];
			contrib = [node insamp model_params' Z'];
			if (j==1) contribs = zeros(n_pooled, columns(contrib)); endif
			contribs(j,:) = contrib;
		endfor
		MPI_Send(contribs, 0, mytag, CW);
		pause(0.1); # give time for the fronted to send a stop message, if done
		# check if we're done
		if (MPI_Iprobe(0, node, CW)) # check for ping from rank 0
			junk = MPI_Recv(0, node, CW);
			break;
		endif
	endwhile
else # frontend
	received = 0;
	done = false;
	while received < reps
		# retrieve results from compute nodes
		pause(0.1);
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
				if verbose
					printf("\nContribution received from node%d.  Received so far: %d\n", i, received);
				endif
				if done
					# tell compute nodes to stop loop
					for j = 1:5
						for i = 1:(nodes-1)
							if (j==1) MPI_Send(" ",i,i,CW); endif # send out message to stop
							ready = MPI_Iprobe(i, mytag, CW); # get last messages
							if ready contribs = MPI_Recv(i, mytag, CW); endif
						endfor
					endfor
					break;
				endif
			endif
		endfor
	endwhile
endif
if not(MPI_Finalized) MPI_Finalize; endif
