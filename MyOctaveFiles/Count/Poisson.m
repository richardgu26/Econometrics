function [logdensity, score, prediction] = Poisson(theta, data)
	y = data(:,1);
	x = data(:,2:columns(data));
	lambda = exp(x*theta);
	logdensity = y .* (x*theta) - lambda - mc_lgamma(y+1);
	e = y - lambda;
	score = diag(e)*x;
#	score = "na"; # uncomment to use numeric score vector
	prediction = lambda;
endfunction
