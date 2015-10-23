# Copyright (C) 2003,2004, 2005  Michael Creel <michael.creel@uab.es>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

# GMM example file, shows initial consistent estimator,
# estimation of efficient weight, and second round
# efficient estimator

1;

# the form a moment condition should take
function m = poisson_moments(theta, data, momentargs)
	k = momentargs{1}; # use this so that data can hold dep, indeps, and instr
	y = data(:,1);
	x = data(:,2:k+1);
	w = data(:, k+2:columns(data));
	lambda = exp(x*theta);
	e = y ./ lambda - 1;
	m = diag(e)*w;
endfunction	


# make regressors
n = 100;
k = 5;
x = [ones(n,1) randn(n,k-1)];
# instruments
w = x;
# generate dep var
theta_true = ones(k,1);
lambda = exp(x*theta_true);
y = poissrnd(lambda);
# The arguments for gmm_estimate
theta = zeros(k,1); # start values
data = [y x w];
weight = eye(columns(w));
moments = "poisson_moments";
momentargs = {k}; # needed to know where x ends and w starts

# initial consistent estimate: only used to get moment covariance (needed for t-stats) no screen output
[theta, obj_value, convergence] = gmm_estimate(theta, data, weight, moments, momentargs);

# moment covariance
# this method is valid when moments are not autocorrelated
# the user is reponsible to properly estimate the efficient weight
m = feval(moments, theta, data, momentargs);
momentcov = cov(m);

# estimation results (no need for efficient weight, exact identification)
gmmtitle = "Poisson GMM example";
gmm_results(theta, data, weight, moments, momentargs);

printf("\nThe true parameter values:\n");
prettyprint_c(theta_true', char(["theta1"; "theta2"; "theta3"; "theta4"; "theta5"]));

