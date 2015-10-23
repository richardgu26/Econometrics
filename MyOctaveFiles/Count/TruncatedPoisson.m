#  Truncated poisson 
function output = TruncatedPoisson(theta, otherargs)
	y = otherargs{1};
	x = otherargs{2};
	k = columns(x);
	lambda = exp(x*theta);
	logdensity =  dmult(y,(x*theta)) - mc_lgamma(y+1) - log(exp(lambda)-1);
	prediction = lambda ./ (1 - exp(-lambda));
	score = y - prediction;
	score = dmult(score,x);
	output = {logdensity, score, prediction};
endfunction
