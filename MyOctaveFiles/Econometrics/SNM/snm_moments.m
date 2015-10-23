function m = snm_moments(params, data, momentargs, doplots)

	if nargin < 4 doplots = false; endif	# set this to true to see kernel fits

	model = momentargs{1}; 		# name of the data generating process (DGP)
	modelargs = momentargs{2}; 	# arguments the DGP needs, in addition to params
	L = momentargs{3};		# number of test functions
	K = momentargs{4};		# number of conditioning variables
	M = momentargs{5};		# number of instruments per endog (same for all)
	do_cond = momentargs{7};
	do_uncond = momentargs{8};

	# sizes of things
	n = rows(data);					# sample size

	# split data into endogs, conditioning variables, and instruments
	instruments = data(:,L+K+1:L+K+M);

	# get residuals (conditional and unconditional)
	[ce ue]  = snm_residuals(params, data, momentargs);

	if do_cond
		# interact instruments with errors
		m = zeros(n, L*M);
		inst_i = instruments;
		for i=1:L
			if (columns(instruments) != M)
				inst_i = instruments(:,i*M-M+1:i*M);
			endif
			m(:,i*M-M+1:i*M) = inst_i .* repmat(ce(:,i), 1, M); # errors interacted with instruments
		endfor
	endif

	if (do_uncond) & (!do_cond)
		m = ue;
	endif

	if do_cond & do_uncond
		m = [m ue];
	endif

endfunction

