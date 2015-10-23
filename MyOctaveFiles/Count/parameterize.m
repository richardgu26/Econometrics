# this is a collection of parameter transformations so that
# estimation routines can apply the Delta Method.

function theta = parameterize(theta, otherargs)
data = otherargs{1};
model = otherargs{2};
modelargs = otherargs{3};
nb_test = strcmp(model, "NegativeBinomial");
nb_test = nb_test + strcmp(model, "TruncatedNegativeBinomial");
nb_test = nb_test + strcmp(model, "NegBin");

# the NB family
if (nb_test > 0)
	p = rows(theta);
	alpha = eps + exp(theta(p,1));
	theta(p, 1) = alpha;
	
# GURMU
elseif strcmp(model, "Gurmu") || strcmp(model, "GurmuBinary") || strcmp(model, "GurmuTruncated")
	# parameters are beta, lnalpha, lnlambda, c.
	k = columns(data) - 1;
	# alpha must be positive, and not too big, so that gamma() doesn't crash
	# this parameterization forces it to stay in bounds
	alpha = 165/(1 + exp(5 - theta(k+1,:))); # eps < alpha < 165
	theta(k + 1, 1) = alpha;

# MNB
elseif strcmp(model, "MixNegativeBinomial")

	# raw parameters are beta, lnalpha, beta2, lnalpha2, mix
	# this makes the 2 alphas > 0, and 0 < mix < 0.5
	k = columns(data) - 1;
  	beta = theta(1:k,:);
  	alpha1 = eps + exp(theta(k+1,:));
  	beta2 = theta(k+2: 2*k+1,:);
	alpha2 = eps + exp(theta(2*k+2,:));
 	mix = theta(2*k+3,:);
	# mix parameter is in (0.01,0.5). For identification, first component is dominant
	# there is a lower bound on mix to prevent loss of identification of second component
	mix = 0.01 + (0.98)/(1 + exp(-1-mix)); # the -2 means that first component starts with high weight
	theta = [beta; alpha1; beta2; alpha2; mix];


# CMNB
elseif strcmp(model, "ConstrainedMixNegativeBinomial")

	# raw parameters are beta, lnalpha, c2, lnalpha2, mix
	# this makes the 2 alphas > 0, and 0 < mix < 0.5
	k = columns(data) - 1;
	const1 = theta(1,:);
  	beta = theta(2:k,:);
  	alpha1 = eps + exp(theta(k+1,:));
	const2 = theta(k+2,:);
	alpha2 = eps + exp(theta(k+3,:));
 	mix = theta(k+4,:);
	# mix parameter is in (0.01,0.5). For identification, first component is dominant
	# there is a lower bound on mix to prevent loss of identification of second component
	mix = 0.01 + (0.98)/(1 + exp(-1-mix)); # the -2 means that first component starts with high weight
	theta = [const1; beta; alpha1; const2; alpha2; mix];


# the C++ version doesn't use this yet
# # NBSNP
# elseif strcmp(model, "NegBinSNP")
# 
# 	# raw parameters are beta, lnalpha, gam_1, .. gam_p
# 	# this makes the alphas > 0
# 	# NegBinSNP internally scales the gammas, this is not
# 	# done here, otherwise, the higher order ones would become
# 	# very close to zeros after scaling
# 	k = columns(data) - 1;
#  	alpha = eps + exp(theta(k+1,:));
# 	theta(k+1,:) = alpha;
# 	p = rows(theta);
# 	for j = k+2:p
# 		theta(j);
# 		theta(j,:) = theta(j,:) / (10^(j-k-1));
# 		theta(j);
# 	endfor 

endif
	
endfunction	 
