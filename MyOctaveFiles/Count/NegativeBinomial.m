#  Notes: the negbin, hurdle negbin and mixed negbin parameterizations
#  follow Deb and Trivedi, J. Appl. Econometrics,
#  V. 12, 1997, pp. 313 - 336.
# 
#  Use negbin_type (3rd element of data list) to chose Type I or II models
# 


function [logdensity, score, prediction, second_raw_moment] = NegativeBinomial(theta, data, otherargs)
	y = data(:,1);
	x = data(:,2:columns(data));
	negbin_type = otherargs{1};
	k = columns(x);
	theta = parameterize(theta, {data, "NegativeBinomial", otherargs});
	alpha = theta(k + 1,1);
	beta = theta(1:k,:);
 	lambda = eps + exp(x*beta);

	if (negbin_type == 1) psi = eps + lambda/alpha;
	else psi = eps + 1/alpha;
	endif

    	logdensity = mc_lgamma(y + psi) - mc_lgamma(psi) - mc_lgamma(y+1) \
				+ psi .* log(psi ./ (lambda + psi)) + y .* log(lambda ./ (lambda + psi));

%  	# now for the derivatives
%  	if (negbin_type == 1)
%  		dpsi_dalpha = - eps -lambda/(alpha^2);
%  		dpsi_dlambda = eps + 1/alpha;
%  	else
%  		dpsi_dalpha = - eps -1/(alpha^2);
%  		dpsi_dlambda = 0;
%  	endif		
%  
%  	dlambda_dbeta = diag(lambda)*x;
%  
%  
%  	dobj_dpsi = digamma(psi + y) \
%  		- digamma(psi) - log(1 + lambda ./ psi) \
%  		+ (lambda ./ (psi + lambda)) \
%  		- (y ./ (lambda + psi));
%  
%  	dobj_dalpha = dobj_dpsi .* dpsi_dalpha * alpha;
%  
%  	dobj_dlambda = y .* psi ./ lambda ./ (lambda + psi) \
%  		- (psi ./ (lambda + psi));
%  
%  	dobj_dbeta = diag(dobj_dpsi .* dpsi_dlambda)*dlambda_dbeta \
%  		+ diag(dobj_dlambda)*dlambda_dbeta;
%  
%  
%   	score = [dobj_dbeta, dobj_dalpha];

 	prediction = lambda;
	
	if (negbin_type == 1)
		var = lambda .*(1 + alpha);
	else 
		var = lambda .* (1 + alpha*lambda);
	endif
	second_raw_moment = var + prediction .^2;
	score = "na";
endfunction
