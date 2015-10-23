## Copyright (C) 2010 Michael Creel
## 
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## SpecTest: modification of MeasurementErrorIV to check the small sample
## distribution of the GMM criterion test n*s(thetahat)
 
## Author: Michael Creel <michael@yosemite>
## Created: 2010-05-05


# this script does a Monte Carlo on the GMM criterion test, using an overidentified IV estimator
# for a simple dynamic model with measurement error.
# Note that this uses the general iterative computation of GMM estimator, thought of as 
# a minimizer. This could be dramatically sped up by using the analytic GIV formula.

1;

function m = moments(theta, data)
	maxlag = 2;
	n = rows(data);
	y = data(maxlag+1:n,1);
	ylag = data(maxlag+1:n,2);
	x = data(maxlag+1:n,3);
	xlags = lags(data(:,3), maxlag);
	xlags = xlags(maxlag+1:n,:);
	n = rows(y);
	regressors = [ones(n,1) ylag x];
	e = y - regressors*theta;
	instruments = [ones(n,1) x xlags];
	m = diag(e)*instruments;
endfunction	
	
# the function to be Monte Carlo'ed
function contrib = wrapper(args)
	n = args{1};    # sample size
	sig = args{2};  # st. dev. of meas. error
	x = randn(n,1); # an exogenous regressor
	e = randn(n,1); # the error term
	ystar = zeros(n,1);
	# generate the dep var
	for t = 2:n
	  ystar(t,:) = 0 + 0.9*ystar(t-1,:) + 1*x(t,:) + e(t,:);
	endfor
    	# add measurement error
	y = ystar + sig*randn(n,1);
	ylag = lag(y,1);
	data = [y ylag x];
	data = data(2:n,:); # drop first obs, missing due to lag
	theta = [0; 0.9; 1];
	weight = 1;

	maxlag = 2;
	n = rows(data);
	y = data(maxlag+1:n,1);
	ylag = data(maxlag+1:n,2);
	x = data(maxlag+1:n,3);
	xlags = lags(data(:,3), maxlag);
	xlags = xlags(maxlag+1:n,:);
	n = rows(y);
	regressors = [ones(n,1) ylag x];
	e = y - regressors*theta;
	instruments = [ones(n,1) x xlags];
	xhat = instruments*ols(regressors,instruments);
	thetahat = inv(regressors'*xhat)*xhat'*y;
	% first round GMM to get initial consistent estimator, for computation of efficient weight matrix
	%thetahat = gmm_estimate(theta, data, weight, "moments");
	% efficient weight matrix
	m = feval("moments", thetahat, data);
	weight = inverse(cov(m));
	# second round efficient GMM
	[thetahat, objvalue] = gmm_estimate(theta, data, weight, "moments");
	objvalue = n*objvalue;
	df = columns(m) - rows(theta);
	reject = objvalue > chi2inv([0.9 0.95 0.99], df); # does the test reject at these signif levels?
	contrib = [objvalue reject];
endfunction

# do the Monte Carlo
n = 50;
sig = 1;
args = {n, sig};
outfile = 'SpecTest.out';
reps = 1000;
system('rm SpecTest.out'); # start from scratch each time
n_pooled = 100;
verbose = true;
montecarlo('wrapper', args, reps, outfile, n_pooled, true);


# analyze results
load SpecTest.out;
results = SpecTest(:,3:6); # drop the timing and node info added by montecarlo.m
clear SpecTest;
close all;

hist(results(:,1),30);
dstats(results);
		


