# Note: for estimation you're better of using
# GurmuBinary+GurmuTruncated: this is faster
# The likelihood fn and caic, etc. are simple the sum of
# the two components. This routine is mostly
# useful for prediction
#
# Note: this fixes restrict = 1. The prediction is NOT valid if you change this
function [logdensity, score, prediction] = GurmuHurdle(theta, data, otherargs)
	y = data(:,1);
	x = data(:,2:columns(data));
	k = columns(x);

	restrict = 1;
	p = rows(theta)/2;
	theta1 = theta(1:p,:);	
	theta1 = parameterize(theta1, {data, "Gurmu", {}});
	theta2 = theta(p+1:2*p,:);
	theta2 = parameterize(theta2, {data, "Gurmu", {}});

# First the Bernoulli part
	test = y == 0;
	prob_0 = GurmuMGF(theta1, x, y-y);
	logdensity1 = test .* log(prob_0) + (1 - test) .* log(1 - prob_0);

# Now the truncated part
	k = columns(x);
	beta = theta2(1:k,:); 
	M1 = GurmuMGF(theta2, x, y);
	M2 = GurmuMGF(theta2, x, y-y);
	prob_pos = 1 - M2;
	logdensity2 =  y.* (x*beta)  - mc_lgamma(y + 1) + log(M1) - log(1 - M2);
  	logdensity2 = logdensity2 .* (1 - test);

	logdensity = logdensity1 + logdensity2;


# prediction	
	lambda = exp(x*beta);
	prediction = lambda .* (1 - prob_0) ./ prob_pos;

	score = "na";
endfunction
