function [endogs, condvars, P] = snm_dataprep(endogs, condvars, P) 
	data = [endogs condvars];
	L = columns(endogs);
	K = columns(condvars);
	if (nargin < 3)
		P = inv(chol(cov([condvars])));
		P = [[diag(1 ./ std(endogs)) zeros(L,K)]; [zeros(K,L) P]]; 
	endif

	data = data*P;
	endogs = data(:,1:L);
	condvars = data(:,L+1:columns(data));
endfunction
