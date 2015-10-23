# Hurdle Negative binomial
# Use negbin_type (3rd element of data list) to chose Type I or II models
function [logdensity, score, prediction] = HurdleNegativeBinomial(theta, data, otherargs)
	y = data(:,1);
	x = data(:,2:columns(data));
	negbin_type = otherargs{1};
	
	k = columns(x);
		
	beta1 = theta(1:k,:);
	beta2 = theta(k+1:2*k,:);
	alpha = eps + exp(theta(2*k+1,:));
	
	# First the Bernoulli part
	lambda = eps + exp(x*beta1);

	if (negbin_type == 1) psi = lambda;  # alpha normalized to one here for identification 
	else psi = 1;
	endif

	test = (y == 0);

	prob_0 = eps + (1-2*eps)*(psi ./ (lambda + psi)) .^ psi;
	
	logdensity1 = test .* log(prob_0) + (1 - test) .* log(1 - prob_0);
	
	# Now the TNB part
	lambda = eps + exp(x*beta2);
	
	if (negbin_type == 1) psi = eps + lambda/alpha;
	else psi = eps + 1/alpha;
	endif
	
	# top-bound psi to avoid crashes
	testpsi = psi > 10000;
	psi = 10000*testpsi + (1 - testpsi).*psi;
	
	
	prob_pos = 1 - (psi ./ (lambda + psi)) .^ psi;
	prob_pos = (1 - 2*eps)*prob_pos + eps;
 
	logdensity2 = mc_lgamma(y + psi) - mc_lgamma(y + 1) - mc_lgamma(psi) \
					+ psi .* log(psi) - psi .* log(lambda + psi)  \
					+ y .* log(lambda) - y .* log(lambda + psi) \
					- log(prob_pos);
	logdensity2 = (1 - test) .* logdensity2;
	
	# Now combine things
	logdensity = logdensity1 + logdensity2;
	prediction = lambda .* (1 - prob_0) ./ prob_pos;
	score = "na";
	
endfunction	
