#  Truncated negative binomial
# Use negbin_type (3rd element of data list) to chose Type I or II models
function [logdensity, score, prediction] = TruncatedNegativeBinomial(theta, otherargs)
	theta = parameterize({theta, "NegativeBinomial", data});

	y = otherargs{1};
	x = otherargs{2};
	negbin_type = otherargs{3};

	k = columns(x);
	beta = theta(1:k,:);
	alpha = theta(k+1,:);

	lambda = exp(x*beta);
	
	if (negbin_type == 1) psi = lambda/alpha;
	else psi = 1/alpha;
	endif
	 
	logdensity = mc_lgamma(y + psi) - mc_lgamma(y + 1) - mc_lgamma(psi) \
					+ psi .* log(psi) - psi .* log(lambda + psi)  \
					+ y .* log(lambda) - y .* log(lambda + psi) \
					- log(1 - (psi ./(lambda + psi)) .^ psi);

	score = "na";
	prob_pos = 1 - ((psi ./ (lambda + psi)) .^ psi);
	prediction = lambda ./ prob_pos;

endfunction
