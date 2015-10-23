% this is the log-likelihood for the linear model
% y = x*b+e, with e~N(0,s^2)
% the parameter theta = [b' s]
% this file is used for demonstration of MLE, to compare to OLS
function logdensity = NormalLF(theta, y, x)
    k = size(theta,1);
    b = theta(1:k-1,:);
    s = theta(k,:);
    e = y - x*b;
    % LF for single each observation
    logdensity = -log(sqrt(2*pi)) - log(s) - e.*e/(2*s^2);
    % sum up
    logdensity = mean(logdensity);
end


