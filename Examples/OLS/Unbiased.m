# this is a little Monte Carlo exercise that illustrates that
# the OLS estimator is unbiased when we have strong exogeneity

reps = 1000; # number of Monte Carlo reps.
n = 20; # sample size
sig = 3;  # st. dev. of errors

x = [ones(n,1), randn(n,1)];  # x is fixed over repeated samples
beta = [1; 2]; # true beta
true = x*beta;

e = sig*randn(n,reps);
y = kron(true, ones(1,reps)) + e;

[betas, sigsq, e] = ols(y,x);

betas(2,:) = betas(2,:) - 2;
title("Beta hat - Beta true");
hist(betas(2,:),30,1);
#legend("off");
pause(10);

%print("Unbiased.svg", "-dsvg");
	
