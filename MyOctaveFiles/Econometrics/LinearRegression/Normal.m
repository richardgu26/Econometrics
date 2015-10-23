% this is the log-likelihood for the linear model
% y = x*b+e, with e~N(0,s^2)
% the parameter theta = [b' s]
% this file is used for demonstration of MLE, to compare to OLS
function [logdensity, score, prediction] = Normal(theta, data)
	[n k] = size(data);
	k = k-1;
	y = data(:,1);
	x = data(:,2:k+1);
	b = theta(1:k,:);
	s = theta(k+1,:);
	e = y - x*b;
	logdensity = -log(sqrt(2*pi)) - log(s) - e.*e/(2*s^2);
	score = 'na'; % to lazy to fill this in!
	prediction = x*b;
end
