#  Mixture Negative binomial (I and II, negbin_type chooses)
#  Deb and Trivedi, J. Appl. Econometrics, V. 12, 1997, pp. 313 - 336.
# 
#  parameters are ordered to use simple NB parameters as start for first
#  component, then get the constant and alpha of second component,
#  and finally the mixture

function [logdensity, score, prediction, second_raw_moment] = ConstrainedMixNegativeBinomial(theta, data, otherargs)
	y = data(:,1);
	x = data(:,2:columns(data));
	negbin_type = otherargs{1};
	n = rows(x);
	k = columns(x);
	
	theta = parameterize(theta, {data, "ConstrainedMixNegativeBinomial", {}});
	
	const1 = theta(1,:);
  	beta = theta(2:k,:);
  	alpha1 = theta(k+1,:);
	const2 = theta(k+2,:);
	alpha2 = theta(k+3,:);
 	mix = theta(k+4,:);

	beta1 = [const1; beta];
	beta2 = [const2; beta];

  	lambda1 = eps + exp(x*beta1);
  	lambda2 = eps + exp(x*beta2);

	if (negbin_type == 1)
		psi1 = eps + lambda1/alpha1;
		psi2 = eps + lambda2/alpha2;
	else 
		psi1 = eps + ones(n,1)/alpha1;
		psi2 = eps + ones(n,1)/alpha2;
	endif

	# top-bound psi to avoid crashes
	testpsi = psi1 > 10000;
	psi1 = 10000*testpsi + (1 - testpsi).*psi1;
	testpsi = psi2 > 10000;
	psi2 = 10000*testpsi + (1 - testpsi).*psi2;


	t1 = lambda1 + psi1;
	t2 = lambda2 + psi2;
	logdensity1 = mc_lgamma(y + psi1) - mc_lgamma(psi1) - mc_lgamma(y+1) \
				+ psi1 .* log(psi1 ./ t1) + y .* log(lambda1 ./ t1);

	logdensity2 = mc_lgamma(y + psi2) - mc_lgamma(psi2) - mc_lgamma(y+1) \
				+ psi2 .* log(psi2 ./ t2) + y .* log(lambda2 ./ t2);

	density1 = exp(logdensity1);
	density2 = exp(logdensity2);

	mix_density = mix*density1 + (1-mix)*density2;
	logdensity = log(eps + mix_density);


# 	# now for the derivatives
# 	if (negbin_type == 1)
# 		dpsi1_dalpha1	 = -lambda1/(alpha1^2);
# 		dpsi1_dlambda1 = 1/alpha1;
# 		dpsi2_dalpha2	 = -lambda2/(alpha2^2);
# 		dpsi2_dlambda2 = 1/alpha2;
# 		
# 	else
# 		dpsi1_dalpha1 = -1/(alpha1^2);
# 		dpsi1_dlambda1 = 0;
# 		dpsi2_dalpha2 = -1/(alpha2^2);
# 		dpsi2_dlambda2 = 0;
# 	endif	
# 	
# 
# 		# first derivatives of log component densities (taken from negbin)
# 		dpsi1 = mc_digamma(psi1 + y) \
# 			- mc_digamma(psi1) - log(1 + lambda1 ./ psi1) \
# 			+ (lambda1 ./ (psi1 + lambda1)) \
# 			- (y ./ (lambda1 + psi1));
# 
# 		dalpha1 = dpsi1 .* dpsi1_dalpha1 * alpha1;
# 
# 		dlambda1 = y .* psi1 ./ lambda1 .* (1 ./ (lambda1 + psi1)) \
# 			- (psi1 ./ (lambda1 + psi1));
# 
# 		dbeta1 = dpsi1 .* dpsi1_dlambda1 .* lambda1 \
# 			+ dlambda1 .* lambda1;
# 		dbeta1 = dmult(dbeta1, x);	
# 
# 		dpsi2 = mc_digamma(psi2 + y) \
# 			- mc_digamma(psi2) - log(1 + lambda2 ./ psi2) \
# 			+ (lambda2 ./ (psi2 + lambda2)) \
# 			- (y ./ (lambda2 + psi2));
# 
# 		dalpha2 = dpsi2 .* dpsi2_dalpha2 * alpha2;
# 
# 		dlambda2 = y .* psi2 ./ lambda2 .* (1 ./ (lambda2 + psi2)) \
# 			- (psi2 ./ (lambda2 + psi2));
# 
# 		dbeta2 = dpsi2 .* dpsi2_dlambda2 .* lambda2 \
# 			+ dlambda2 .* lambda2;
# 		dbeta2 = dmult(dbeta2, x);	
# 			
#   	beta = theta(1:k,:);
			
# 		# this is the relationship between log component and log mixture densities
# 		dbeta1 = dmult(mix .* density1 ./ mix_density, dbeta1);
# 		dalpha1 = dalpha1 ./ mix_density .* mix .* density1;
# 		dbeta2 = dmult(mix .* density2 ./ mix_density, dbeta2);
# 		dalpha2 = dalpha2 ./ mix_density .* (1-mix) .* density2;
# 	
# 		dbeta = dbeta1(:,2:k) + dbeta2(:,2:k);
# 		
# 		dconst1 = dbeta1(:,1);
# 		dconst2 = dbeta2(:,1);
# 		
# 		dmix = (density1 - density2) ./ mix_density;
# 		dmix = dmix .* (mix .^ 2) * exp(-theta(k+4,:)); # adjust to get deriv wrt raw param
# 		
# 		score = [dconst1 dbeta , dalpha1 , dconst2 , dalpha2 , dmix];
		score = "na";
		prediction = mix * lambda1 + (1 - mix) * lambda2;
		if (negbin_type == 1)
			var1 = lambda1 .* (1 + alpha1);
			var2 = lambda2 .* (1 + alpha2);
		else 
			var1 = lambda1 .* (1 + alpha1*lambda1);
			var2 = lambda2 .* (1 + alpha2*lambda2);
		endif
		var = mix*var1 + (1 - mix) * var2;
		second_raw_moment = var + prediction .^2;
endfunction
