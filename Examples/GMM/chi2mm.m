# method of moments estimators for a sample from Chi^2(theta)
# increase n to observe consistency
# note that the two estimators are different from one another
n = 30;
theta = 3;
y=chi2rnd(theta, n, 1);
thetahat = mean(y);
thetahat

thetahat = 0.5*var(y);
thetahat

