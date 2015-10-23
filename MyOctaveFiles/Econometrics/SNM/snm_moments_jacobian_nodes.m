function D = snm_moments_jacobian_nodes(args);

	theta = args{1};
	n = args{2};
	momentargs = args{3};

	model = momentargs{1};
	modelargs = momentargs{2};

	# get bootstrap data (size n)
	modelargs{1} = "";
	modelargs{2} = n;
	modelargs{4} = true; # make instruments
	[data, L, K, snm_draws] = feval(model, theta, modelargs);
	M = columns(data) - L - K;

	modelargs{1} = snm_draws;

	# split data into endogs, conditioning variables, and instruments
	endogs = data(:,1:L);
	condvars = data(:,L+1:L+K);
	[endogs, condvars, P] = snm_dataprep(endogs, condvars);
	instruments = data(:,L+K+1:L+K+M);
	data = [endogs, condvars, instruments];

	D = numgradient("average_moments", {theta, data, "snm_moments", momentargs, 0});
	D = D(:)';

endfunction

