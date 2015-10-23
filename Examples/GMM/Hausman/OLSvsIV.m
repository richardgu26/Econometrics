# this is a little Monte Carlo exercise that illustrates that
# the OLS estimator is biased and inconsistent when errors are
# correlated with regressors, but that the IV estimator is consistent

reps = 1000; # number of Monte Carlo reps.
betaols = zeros(reps,2);
betaiv = betaols;
n = 1000; # sample size

# covariance of X, W, e
cov_X_W = 0.2;  % experiment with lowering or raising this: quality of instrument
cov_X_e = 0;
sig = [
3, cov_X_W, cov_X_e;
cov_X_W, 1, 0;
cov_X_e, 0, 1];

true = [1; 2]; # true beta
p = chol(sig);
for i = 1:reps
	x = randn(n,3)*p;
	e = x(:,3);
	w = [ones(n,1), x(:,2)];
	x = [ones(n,1), x(:,1)];
	y = x*true + e;
	# OLS
	[b, sigsq, e] = ols(y,x);
	betaols(i,:) = betaols(i,:) + b';
	# IV
	b = inv(w'*x)*w'*y;
	betaiv(i,:) = betaiv(i,:) + b';

endfor
betaols = betaols';
betaiv = betaiv';
close;

hist(betaols(2,:), 50, 1);
title("OLS estimates");
%print("ols.svg", "-dsvg");
figure;
hist(betaiv(2,:), 50, 1);
title("IV estimates");
%print("iv.svg", "-dsvg");

printf("OLS results\n");
dstats(betaols');
printf("IV results\n");
dstats(betaiv');
