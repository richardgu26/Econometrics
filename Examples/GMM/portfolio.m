# 
# This is an example of estimation of a portfolio model
# as in Hansen and Singleton, 1982. The data comes from
# Tauchen, JBES, 1986 (available by ftp).
# The notation follows my lecture notes. In particular
# gamma is the coefficient of RRA. Try out different sets
# of intruments and different lag lengths. You will see
# that beta is stable across  models, but that gamma is not.
# 
# Michael Creel, Oct. 24, 2000
# revised 7 Nov. 2002
# translated to Octave 9/9/2003
# michael.creel@uab.es


# The main thing you have to to is define the moments for estimation,
# as in the function that follows immediately, and ensure that the 
# efficient weight matrix is estimated appropriately. The rest
# is already programmed.


# the following function defines the moment conditions, given the parameter
# vector and the instruments. It returns the individual contributions
# to allow estimation of the covariance matrix of the moments.
1;
function m = portfolio_moments(theta, data, momentargs)

	# parameters
	beta = theta(1,1);
	gam = theta(2,1);
	
	# data items
	c = data(:,1);
	r = data(:,2);
	inst = data(:,3:columns(data));
	
	#  form error function
	# note that c = c_t / c_t-1 (for stationarity) was done in data preparation
	e = 1 - beta*(1 + r) .* (c .^ (-gam)); 

	# cross with instruments
	m = diag(e)*inst;
endfunction

	
data = load("tauchen.data");

c = data(:,1);
p = data(:,2);
d = data(:,3);


# form net return and stationary consumption
r = (p + d) ./ lag(p,1) - 1;
c = c ./ lag(c,1); # ensure stationarity

# choose maximal lag of instruments.
max_lag = 1;
inst = [c r d p];
inst = lags(inst, max_lag);
inst = st_norm(inst);
inst = [ones(rows(inst),1) inst];
# drop rows with missing values
data = [c r inst];
data = data(max_lag+1:rows(data),:);

# now make things the way the moment conditions function likes them
momentargs = {};
theta = zeros(2,1);
moments = "portfolio_moments";	
weight = eye(columns(inst));
names = char("beta","gamma");
title = "Example of GMM estimation of rational expectations model";

# initial consistent estimate
theta = [0.9;1];
ub = [1;5];
lb = [0.8; -5];
nt = 3;
ns = 1;
rt = 0.5; # careful - this is too low for many problems
maxevals = 1000;
neps = 5;
functol = 1e-10;
paramtol = 1e-3;
verbosity = 2; # only final results. Inc
minarg = 1;
control = { lb, ub, nt, ns, rt, maxevals, neps, functol, paramtol, verbosity, 1};
% some SA to get good start values
[theta, obj_value, convergence] = gmm_estimate(theta, data, weight, moments, momentargs, control);
% BFGS to finish
[theta, obj_value, convergence] = gmm_estimate(theta, data, weight, moments, momentargs);
theta

# efficient weight matrix
m = feval(moments, theta, data, momentargs);
momentcov = cov(m);
weight = inverse(momentcov);
% some SA to get good start values
[theta, obj_value, convergence] = gmm_estimate(theta, data, weight, moments, momentargs, control);
% BFGS finish, and report results
gmm_results(theta, data, weight, moments, momentargs, names, title);





