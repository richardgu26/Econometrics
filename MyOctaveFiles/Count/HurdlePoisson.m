#  Hurdle Poisson 
#  All in one, rather than binary followed by truncated Poisson
# (not computationally efficient, but easier to deal with)
function [logdensity, score, prediction] = HurdlePoisson(theta, otherargs)
	y = otherargs{1};
	x = otherargs{2};
	k = columns(x);
	test = y == 0;
	
	# First the Bernoulli part
	beta = theta(1:k,:);
	lambda = exp(x*beta);
	prob_0 = exp(-lambda);
	logdensity1 = test .* log(prob_0) + (1-test) .* log(1 - prob_0);
	
	# Now the truncated part
	beta = theta(k+1:rows(theta),:);
	lambda = exp(x*beta);
	prob_pos = 1 - exp(-lambda);
	logdensity2 = dmult(y,(x*beta)) - mc_lgamma(y+1) - log(exp(lambda)-1); 
	logdensity2 = (1 - test) .* logdensity2;
	# Conbine things
	logdensity = logdensity1 + logdensity2;
	prediction = lambda .* (1 - prob_0) ./ prob_pos;
	score = "na";
endfunction
