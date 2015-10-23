# Copyright (C) 2003,2004,2005  Michael Creel <michael.creel@uab.es>
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

# usage: [theta, V, obj_value] =
#  gmm_results(theta, data, weight, moments, momentargs, names, title, unscale, control, compute_nodes, momentcov)
#
# inputs:
#      theta: column vector initial parameters
#       data: data matrix
#     weight: the GMM weight matrix
#    moments: name of function computes the moments
#             (should return nXg matrix of contributions)
# momentargs: (cell) additional inputs needed to compute moments.
#             May be empty ("")
#      names: vector of parameter names
#             e.g., names = char("param1", "param2");
#      title: string, describes model estimated
#    unscale: (optional) cell that holds means and std. dev. of data
#             (see scale_data)
#    control: (optional) BFGS or SA controls (see bfgsmin and samin). May be empty ("").
#    compute_nodes: (optional) number of compute nodes if executed in parallel
#             (requires MPITB)
#  momentcov: (optional) Estimated covariance of moment conditions, scaled to converge to finite matrix
#             If this is not supplied, it is assumed that the weight matrix is efficient,
# 	      so that momentcov = inv(weight)
# outputs:
# theta: GMM estimated parameters
# V: estimate of covariance of parameters
# obj_value: the value of the GMM objective function
#
# please type "gmm_example" while in octave to see an example


function [theta, V, obj_value] = gmm_results(theta, data, weight, moments, momentargs, names, title, unscale, control, compute_nodes, momentcov)

	efficient = false;
	if (nargin < 11) efficient = true; endif
	if (nargin == 11) && (weight == inv(momentcov)) efficient = true; endif
	if (nargin < 10) compute_nodes = 0; endif # serial by default
	if (nargin < 9) control = ""; endif
	if (nargin < 7) title = "GMM results"; endif
	if (nargin < 6)
		names = 1:rows(theta);
		names = names';
	end

	[theta, obj_value, convergence] = gmm_estimate(theta, data, weight, moments, momentargs, control, compute_nodes);


	m = feval(moments, theta, data, momentargs); # find out how many obsns. and moments we have
	n = rows(m);
	q = columns(m);

	if convergence == 1
		convergence="Normal convergence";
	else
		convergence="No convergence";
	endif

	if efficient
		V = gmm_variance(theta, data, weight, moments, momentargs);
	else
		V = gmm_variance_inefficient(theta, data, weight, momentcov, moments, momentargs);
	endif
	# unscale results if argument has been passed
	# this puts coefficients into scale corresponding to the original data
	if nargin > 7
		if iscell(unscale)
			[theta, V] = unscale_parameters(theta, V, unscale);
		endif
	endif

	[theta, V] = delta_method("parameterize", theta, {data, moments, momentargs}, V);

	k = rows(theta);
	se = sqrt(diag(V));

	printf("\n\n******************************************************\n");
	disp(title);
	printf("\nGMM Estimation Results\n");
	printf("BFGS convergence: %s\n", convergence);
	printf("\nObjective function value: %f\n", obj_value);
	printf("Observations: %d\n", n);
	if (efficient && (nargin == 11)) printf("Using efficient weight matrix\n"); endif
	if (efficient && (nargin < 11)) printf("No moment covariance supplied, assuming efficient weight matrix\n"); endif
	junk = "X^2 test";
	df = q - k;
	if ((df > 0) && efficient)
		clabels = char("Value","df","p-value");
		a = [n*obj_value, df, 1 - chi2cdf(n*obj_value, df)];
		printf("\n");
		prettyprint(a, junk, clabels);
	endif
	if (df == 0)
		disp("\nExactly identified, no spec. test");
	end;

	if !efficient
		disp("\nUsing inefficient weight, no spec. test");
	endif

	# results for parameters
	a =[theta, se, theta./se, 2 - 2*normcdf(abs(theta ./ se))];
	clabels = char("estimate", "st. err", "t-stat", "p-value");
	printf("\n");
	prettyprint(a, names, clabels);

	printf("******************************************************\n");
endfunction
