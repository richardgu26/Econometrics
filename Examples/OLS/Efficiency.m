# this is a little Monte Carlo exercise that illustrates that
# the OLS estimator is efficient. It compares the OLS estimator
# to an estimator that is the average of the OLS estimators using 
# 3 subsamples

reps = 10000; # number of Monte Carlo reps.
n = 21; # sample size
sig = 2;  # st. dev. of errors

x = [ones(n,1), randn(n,1)];  # x is fixed over repeated samples
beta = [1; 2]; # true beta
true = x*beta;

betas = zeros(reps,2); # holder for results

for i = 1:reps

	# generate dep var
	e = sig*randn(n,1);
	y = true + e;

	# the OLS estimator
	[beta_ols, sigsq, e] = ols(y,x);
	betas(i,1) = beta_ols(2,1);

	# the average of split sample estimator
	y1 = y(1:7,:);
	y2 = y(8:14,:);
	y3 = y(15:21,:);
	x1 = x(1:7,:);
	x2 = x(8:14,:);
	x3 = x(15:21,:);
	[beta_ss1, sigsq, e] = ols(y1,x1);
	[beta_ss2, sigsq, e] = ols(y2,x2);
	[beta_ss3, sigsq, e] = ols(y3,x3);
	beta_ss = (beta_ss1 + beta_ss3 + beta_ss3)/3; # average the 3 estimators
	betas(i,2) = beta_ss(2,1);

endfor	

hist(betas(:,1),30,1);
legend("off");
title("Beta 2 hat, OLS");
axis([0,4]);
print("efficiency-1.png","-dpng");

hist(betas(:,2),30,1);
title("Beta 2 hat, Split Sample Estimator");
legend("off");
axis([0,4]);
print("efficiency-2.png","-dpng");

	
