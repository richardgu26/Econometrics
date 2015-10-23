# Gurmu semiparametric count model, standard version
# Note that restrict = 1. Prediction NOT valid if you change this
# The standard Gurmu semiparametric model (likelihood function between Gurmu's eqns 9 and 10)
function [logdensity, score, prediction] = Gurmu(theta, data, otherargs)
	y = data(:,1);
	x = data(:,2:columns(data));
	k = columns(x);
	theta = parameterize(theta, {data, "Gurmu", {}});
	beta = theta(1:k,:); 
	logdensity = y .* (x*beta) - mc_lgamma(y + 1) + log(GurmuMGF(theta, x, y));
	score = "na";
	prediction = exp(x*beta); # mean, using restriction
endfunction 
