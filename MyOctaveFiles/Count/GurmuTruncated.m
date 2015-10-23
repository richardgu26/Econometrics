# Note  restrict = 1. Prediction is NOT valid if you change this
function [logdensity, score, prediction] = GurmuTruncated(theta, data, otherargs)
	y = data(:,1);
	x = data(:,2:columns(data));
	k = columns(x);
	theta = parameterize(theta, {data, "Gurmu", {}});
	restrict = 1;
	test = y == 0;
	k = columns(x);
	beta = theta(1:k,:); 
	M1 = GurmuMGF(theta, x, y);
	M2 = GurmuMGF(theta, x, y-y);
	prob_pos = 1 - M2;
	logdensity =  y .* (x*beta)  - mc_lgamma(y + 1) + log(M1) - log(1 - M2);
	logdensity = logdensity .* (1 - test);
	

# prediction	
	lambda = exp(x*beta);
	prediction = lambda ./ prob_pos;

	score = "na";
endfunction
