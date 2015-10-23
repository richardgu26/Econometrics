function fit = jackknife(y, x)
#  this gives jackknife (leave out own observation) 
#  linear fit for y conditional on x
#  useful to check for influential observations
	n = rows(y);
	P = x*inv(x'*x)*x';  # this can use a lot of memory
	M = eye(n) - P;
	h = getdiag(P);
	ols_fit = P*y;
	e = M*y;
	fit = ols_fit - h ./ (1-h) .* e;
endfunction
