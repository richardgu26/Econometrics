# Note restrict = 1. Prediction NOT valid if you change this
function [logdensity, score, prediction] = GurmuBinary(theta, data, otherargs)
	y = data(:,1);
	x = data(:,2:columns(data));
	k = columns(x);
	theta = parameterize(theta, {data, "Gurmu", {}});
	restrict = 1;

# First the Bernoulli part
	test = y == 0;
	prob_0 = GurmuMGF(theta, x, y-y);
	logdensity = test .* log(prob_0) + (1 - test) .* log(1 - prob_0);

# prediction	
	prediction = (1 - prob_0);
	score = "na";
endfunction
