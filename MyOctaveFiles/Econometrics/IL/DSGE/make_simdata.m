% two shock model for IL paper
% NOTE: GetNewBounds is commented out! See at bottom


% there are 2 designs
design = 1;  % set to 1 or 2

n_pooled = 1;   % number of runs accumulated before sent from nodes to frontend
verbose = 1;


if not(MPI_Initialized) MPI_Init; endif
CW = MPI_Comm_Load("NEWORLD");
node = MPI_Comm_rank(CW);
nodes = MPI_Comm_size(CW);
mytag = 48;

if design ==1
 	% design 1
	lb_param_ub = [
	0.20   	0.33   	0.4		% alpha
	0.90    0.99   	0.999		% beta
	0.005    0.025   0.1		% delta
	0     	2      	3		% gam
	0.0    	0.9   	0.999		% rho1
	0.001   0.01 	0.1		% sigma1
	0.0    	0.7     0.999		% rho2
	0.001	0.005  	0.1		% sigma2
	6/24    8/24	12/24		% nss
	];
else design==2
	% design 2
	lb_param_ub = [
	0.20   	0.25   	0.4		% alpha
	0.90    0.97   	0.999		% beta
	0.005    0.04   0.1		% delta
	0     	1      	3		% gam
	0.0    	0.8   	0.999		% rho1
	0.001   0.02 	0.1		% sigma1
	0.0    	0.6     0.999		% rho2
	0.001	0.01  	0.1		% sigma2
	6/24    8/24	12/24		% nss
	];
endif

model_params0 = lb_param_ub(:,2);
lb = lb_param_ub(:,1);
ub = lb_param_ub(:,3);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% you do not need to alter anything below this line %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

first_time = true;

if node
	more_please = 1;
	while more_please
		% break it up to do intermediate writes
		for j = 1:n_pooled
			% re-establish the state of the RNG, dynare seems to synchronize this across MPI nodes
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

			first_time = false; % from now on we want to keep the RNG sequence moving along properly

			% most draws will be from param space, but some from true values
			if outsamp
				model_params = model_params0;
			else
				model_params = rand(rows(ub),1).*(ub-lb) + lb;
			endif

			% break into pieces
			alpha = model_params(1,:);
			beta  = model_params(2,:);
			delta = model_params(3,:);
			gam   = model_params(4,:);
			rho1   = model_params(5,:);
			sigma1 = model_params(6,:);
			rho2   = model_params(7,:);
			sigma2 = model_params(8,:);
			nss   = model_params(9,:);

			% the psi implied by other parameters
			c1 = ((1/beta + delta - 1)/alpha)^(1/(1-alpha));
			kss = nss/c1;
			css = kss * (c1^(1-alpha) - delta);
			c2 = (css)^(-gam/alpha);
			psi = (1-alpha)*((c1/c2)^(-alpha));

			% there are a number of .mod files that each read a different parameter file
			% this writes out the different parameter files so that a number of solutions
			% can be done in parallel
			s = rand("state");
			if     node==1
				save parameterfile1  alpha beta delta gam rho1 sigma1 rho2 sigma2 nss;
				save state1 s;
			elseif node==2
				save parameterfile2   alpha beta delta gam rho1 sigma1 rho2 sigma2 nss;
				save state2 s;
			elseif node==3
				save parameterfile3   alpha beta delta gam rho1 sigma1 rho2 sigma2 nss;
				save state3 s;
			elseif node==4
				save parameterfile4   alpha beta delta gam rho1 sigma1 rho2 sigma2 nss;
				save state4 s;
			elseif node==5
				save parameterfile5   alpha beta delta gam rho1 sigma1 rho2 sigma2 nss;
				save state5 s;
			elseif node==6
				save parameterfile6   alpha beta delta gam rho1 sigma1 rho2 sigma2 nss;
				save state6 s;
			elseif node==7
				save parameterfile7   alpha beta delta gam rho1 sigma1 rho2 sigma2 nss;
				save state7 s;
			elseif node==8
				save parameterfile8   alpha beta delta gam rho1 sigma1 rho2 sigma2 nss;
				save state8 s;
			elseif node==9
				save parameterfile9   alpha beta delta gam rho1 sigma1 rho2 sigma2 nss;
				save state9 s;
			elseif node==10
				save parameterfile10  alpha beta delta gam rho1 sigma1 rho2 sigma2 nss;
				save state10 s;
			elseif node==11
				save parameterfile11  alpha beta delta gam rho1 sigma1 rho2 sigma2 nss;
				save state11 s;
			elseif node==12
				save parameterfile12  alpha beta delta gam rho1 sigma1 rho2 sigma2 nss;
				save state12 s;
			elseif node==13
				save parameterfile13  alpha beta delta gam rho1 sigma1 rho2 sigma2 nss;
				save state13 s;
			elseif node==14
				save parameterfile14  alpha beta delta gam rho1 sigma1 rho2 sigma2 nss;
				save state14 s;
			elseif node==15
				save parameterfile15  alpha beta delta gam rho1 sigma1 rho2 sigma2 nss;
				save state15 s;
			elseif node==16
				save parameterfile16  alpha beta delta gam rho1 sigma1 rho2 sigma2 nss;
				save state16 s;
			elseif node==17
				save parameterfile17  alpha beta delta gam rho1 sigma1 rho2 sigma2 nss;
				save state17 s;
			elseif node==18
				save parameterfile18  alpha beta delta gam rho1 sigma1 rho2 sigma2 nss;
				save state18 s;
			elseif node==19
				save parameterfile19  alpha beta delta gam rho1 sigma1 rho2 sigma2 nss;
				save state19 s;
			elseif node==20
				save parameterfile20  alpha beta delta gam rho1 sigma1 rho2 sigma2 nss;
				save state20 s;
			elseif node==21
				save parameterfile21  alpha beta delta gam rho1 sigma1 rho2 sigma2 nss;
				save state21 s;
			elseif node==22
				save parameterfile22  alpha beta delta gam rho1 sigma1 rho2 sigma2 nss;
				save state22 s;
			elseif node==23
				save parameterfile23  alpha beta delta gam rho1 sigma1 rho2 sigma2 nss;
				save state23 s;
			elseif node==24
				save parameterfile24  alpha beta delta gam rho1 sigma1 rho2 sigma2 nss;
				save state24 s;
			elseif node==25
				save parameterfile25  alpha beta delta gam rho1 sigma1 rho2 sigma2 nss;
				save state25 s;
			elseif node==26
				save parameterfile26  alpha beta delta gam rho1 sigma1 rho2 sigma2 nss;
				save state26 s;
			elseif node==27
				save parameterfile27  alpha beta delta gam rho1 sigma1 rho2 sigma2 nss;
				save state27 s;
			elseif node==28
				save parameterfile28  alpha beta delta gam rho1 sigma1 rho2 sigma2 nss;
				save state28 s;
			elseif node==29
				save parameterfile29  alpha beta delta gam rho1 sigma1 rho2 sigma2 nss;
				save state29 s;
			elseif node==30
				save parameterfile30  alpha beta delta gam rho1 sigma1 rho2 sigma2 nss;
				save state30 s;
			elseif node==31
				save parameterfile31  alpha beta delta gam rho1 sigma1 rho2 sigma2 nss;
				save state31 s;
			elseif node==32
				save parameterfile32  alpha beta delta gam rho1 sigma1 rho2 sigma2 nss;
				save state32 s;
			endif

			% solve model, using the appropriate parameter file
			if node==1
				dynare SimpleModel1 noclearall;
			elseif node==2
				dynare SimpleModel2 noclearall;
			elseif node==3
				dynare SimpleModel3 noclearall;
			elseif node==4
				dynare SimpleModel4 noclearall;
			elseif node==5
				dynare SimpleModel5 noclearall;
			elseif node==6
				dynare SimpleModel6 noclearall;
			elseif node==7
				dynare SimpleModel7 noclearall;
			elseif node==8
				dynare SimpleModel8 noclearall;
			elseif node==9
				dynare SimpleModel9 noclearall;
			elseif node==10
				dynare SimpleModel10 noclearall;
			elseif node==11
				dynare SimpleModel11 noclearall;
			elseif node==12
				dynare SimpleModel12 noclearall;
			elseif node==13
				dynare SimpleModel13 noclearall;
			elseif node==14
				dynare SimpleModel14 noclearall;
			elseif node==15
				dynare SimpleModel15 noclearall;
			elseif node==16
				dynare SimpleModel16 noclearall;
			elseif node==17
				dynare SimpleModel17 noclearall;
			elseif node==18
				dynare SimpleModel18 noclearall;
			elseif node==19
				dynare SimpleModel19 noclearall;
			elseif node==20
				dynare SimpleModel20 noclearall;
			elseif node==21
				dynare SimpleModel21 noclearall;
			elseif node==22
				dynare SimpleModel22 noclearall;
			elseif node==23
				dynare SimpleModel23 noclearall;
			elseif node==24
				dynare SimpleModel24 noclearall;
			elseif node==25
				dynare SimpleModel25 noclearall;
			elseif node==26
				dynare SimpleModel26 noclearall;
			elseif node==27
				dynare SimpleModel27 noclearall;
			elseif node==28
				dynare SimpleModel28 noclearall;
			elseif node==29
				dynare SimpleModel29 noclearall;
			elseif node==30
				dynare SimpleModel30 noclearall;
			elseif node==31
				dynare SimpleModel31 noclearall;
			elseif node==32
				dynare SimpleModel32 noclearall;
			endif
			% get a simulation of length 160 (40 years quarterly), and compute aux. statistic
			data = [y c n MPK MPL];
			data = data(101:260,:);
			Z = aux_stat(data);
			contrib = [outsamp model_params' psi Z'];
			if (i==1) && (j==1) contribs = zeros(n_pooled, columns(contrib)); endif
			contribs(j,:) = contrib;
		endfor
		MPI_Send(contribs, 0, mytag, CW);
		% check if we're done
		if (MPI_Iprobe(0, mytag+1, CW))
			junk = MPI_Recv(0, mytag+1, CW);
			break;
		endif
	endwhile
else % frontend
	received = 0;
	rescaled = false;
	done = false;
	if (exist(sprintf("./%s",outfile)) !=2) recordsize = -1; else recordsize = 1; endif # if output file doesn't exist, we don't know the size of contribs yet
	while received < reps
		% retrieve results from compute nodes
		pause(0.01);
		for i = 1:nodes-1
			% compute nodes have results yet?
			ready = false;
			ready = MPI_Iprobe(i, mytag, CW); % check if message pending
			if ready
				% get it if it's there
				contribs = MPI_Recv(i, mytag, CW);
				% convert to single
				contribs = single(contribs);
				% the first time, write the record size as first entry in data file
				if (recordsize==-1)
					FN = fopen (outfile, "ab");
					if (FN < 0) error ("make_simdata: couldn't open output file %s", outfile); endif
					recordsize = columns(contribs);
					fwrite(FN, recordsize, 'single');
					fclose(FN);
				endif
				need = reps - received;
				received = received + n_pooled;
				% truncate?
				if n_pooled  >= need
						contribs = contribs(1:need,:);
						done = true;
				endif
				% write the data to file
				contribs = contribs';
				FN = fopen (outfile, "ab");
				fwrite(FN, contribs, 'single')
				fclose(FN);
				if verbose
					printf("\nContribution received from node%d.  Received so far: %d\n", i, received);
				endif
				if done
					% tell compute nodes to stop loop
					for j = 1:5
						for i = 1:(nodes-1)
							if (j==1) MPI_Send(" ",i, mytag+1,CW); endif % send out message to stop
							ready = MPI_Iprobe(i, mytag, CW); % get last messages
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
