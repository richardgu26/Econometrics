# Estimates the basic Nerlove Cobb-Douglas model by MLE
# this is intended to show that MLE with normality is the
# ML estimator (for the b's, not for sig^2)

load nerlove.data;

data = data(:,2:6);
data = log(data);
n = rows(data);
y = data(:,1);
x = data(:,2:5);
x = [ones(n,1), x];

# prepare the inputs for mle_results - the first are the starting value for the parameters
# here are some good start values - those from OLS
[theta junk junk ess] = mc_ols(y,x, "", true);
theta = [theta; sqrt(ess/(n-5))];

# if you like, try some other start values, and notice the potential problems!
# theta = [zeros(5,1); 0.5]; 	# initial values, 0 for b's, 1 for s


model = "Normal";			# name of function for log-likelihood
data = [y x]; 				# the data
modelargs = ""; 			# nothing needed here
names = char('constant', 'output', 'labor', 'fuel', 'capital', 'sig'); # parameter names
title = "check MLE with normality, compare to OLS";
mle_results(theta, data, model, modelargs, names, title);
