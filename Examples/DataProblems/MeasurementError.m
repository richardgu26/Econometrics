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

## MeasurementError

## Author: Michael Creel <michael@yosemite>
## Created: 2010-03-04


# this script does a Monte Carlo that shows the inconsistency of OLS estimator
# when there is measurement error of the regressors

1;
# the function to be Monte Carlo'ed
function b = wrapper(args)
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
	# now do OLS, using the data with meas. error in both dep. var. and regressor
	ylag = lag(y,1);
	data = [y ylag x];
	data = data(2:n,:); # drop first obs, missing due to lag
	y = data(:,1);
	ylag = data(:,2);
	x = data(:,3);
	x = [ones(rows(x),1) ylag x]; # data matrix for ols
	b = ols(y, x)';
	b = b - [0 0.9 1]; # subtract true values, so mean should be approx. zero if consistent
endfunction

# do the Monte Carlo
n = 100;
sig = 1;
args = {n, sig};
outfile = 'meas_error.out';
reps = 1000;
system('rm meas_error.out'); # start from scratch each time
n_pooled = 100;
verbose = true;
montecarlo('wrapper', args, reps, outfile, n_pooled, false, verbose);


# analyze results
load meas_error.out;
results = meas_error(:,3:5); # drop the timing and node info added by montecarlo.m
close all;

hist(results(:,1),30);
# print('constant_n100.png', '-dpng');
figure;
hist(results(:,2),30);
%print('ylag_n100.png', '-dpng');
figure
hist(results(:,3),30);
# print('x_n100.png', '-dpng');

dstats(results);
		


