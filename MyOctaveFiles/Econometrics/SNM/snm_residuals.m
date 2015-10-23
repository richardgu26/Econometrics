# ce are conditional errors
# ue are unconditional errors
function [ce, ue] = snm_residuals(params, data, momentargs)

	if (nargin < 7) makevec = false; endif # allow vec-ing out when needed for covariance estimation

	model = momentargs{1}; 		# name of the data generating process (DGP)
	modelargs = momentargs{2}; 	# arguments the DGP needs, in addition to params
	P = momentargs{6};
	do_cond = momentargs{7};
	do_uncond = momentargs{8};


	# get simulated data
	[simdata, L, K] = feval(model, params, modelargs);
	M = columns(simdata) - L - K;

	k = rows(params);
	wwinv = params(k-K+1:k,1); # last L are wws


	# real data
	endogs = data(:,1:L);
	condvars = data(:,L+1:L+K);

	# sizes of things
	n = rows(endogs);  # sample size
	S = rows(simdata);

	# split simulated data into endogs and conditioning variables, and instruments in case of real data
	endogs_sim = simdata(:,1:L);
	condvars_sim = simdata(:,L+1:L+K);
	[endogs_sim, condvars_sim] = snm_dataprep(endogs_sim, condvars_sim, P);
	
	condvars = condvars*diag(wwinv);
	condvars_sim = condvars_sim*diag(wwinv);

	if do_cond
		W = __kernel_weights(condvars_sim, condvars, "__kernel_epanechnikov");
		den = sum(W,2)+eps; # avoid div by zero
		W = diag(1./den)*W;
		# all fits at once - pretty neat!
		fits = W*endogs_sim;
		ce = endogs - fits;
	else ce = "na";
	endif
	if do_uncond
		ue = endogs - repmat(mean(endogs_sim), n, 1);
	else ue = "na";
	endif
endfunction

