# DGP generates probit data
function y = ProbitDGP(theta, x, e, smooth)
	if nargin == 2
		prob1 = normcdf(x*theta);
		n = rows(x);
		y = rand(n,1) < prob1;
	else
		if nargin < 4 error("need 4 args to generate smooth data"); endif
		y = normcdf(smooth*(x*theta - e));
	endif
endfunction
