# This file illustrates how the delta method may be used
# to calculate the covariance of a nonlinear function of
# estimated parameters
1;
# this function that gives the elasticities of x*b wrt x:
function elasticities = ElasticityLinearModel(theta, otherargs)
	x = otherargs{1};
	elasticities = (theta .* x) / (x'*theta);
endfunction


# define some random data	
n = 100;
k = 5;
x = [ones(n,1) rand(n,k-1)]; 
beta = (1:k)'; # true betas are [1,2,3,4,5]'
y = x*beta + randn(n,1);
[b, varb] = mc_ols(y,x);

# we'll evaluate elasticities at the sample means of the data
x = mean(x)';

# the function DeltaMethod numerically differentiates the above function to get the cov. matrix
[elasticities, var_elasticities] = delta_method("ElasticityLinearModel", b, {x}, varb);

# Taa daa - the results!
printf("\n");
printf("Delta method elasticities and standard errors,\nat the means of the regressors\n");
disp([elasticities, diag(var_elasticities)]);
