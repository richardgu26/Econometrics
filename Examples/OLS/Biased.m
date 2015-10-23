% this is a little Monte Carlo exercise that illustrates that
% the OLS estimator is biased when we have an autoregressive model
% with weak exogeneity

reps = 1000; % number of Monte Carlo reps.
n = 20; % sample size
x0 = 0;
betas = zeros(reps,1);
truebetas = [0 0.9];
for i = 1:reps
	x = zeros(n+1,1);
	x(1,1) = 0;

	% generate AR(1) data
	for t = 2:n+1;
		x(t,1) = truebetas(:,1) + truebetas(:,2)*x(t-1) + randn(1,1);
		end
	y = x(2:n+1,1);    % dependent variable
	x = x(1:n,1);      % explanatory variable is the lagged dep var.   
	x = [ones(n,1) x];
	beta = regress(y,x);
	betas(i,1) = beta(2,1);
	end	

betas = betas - truebetas(1,2);
title('Beta hat - Beta true');
hist(betas,30,1);
legend('off');
print('Biased.svg', '-dsvg');
	
