# Copyright (C) 2008 Michael Creel <michael.creel@uab.es>
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

#
# snm_monents_jacobian: estimates Jacobian of average moments (dm/dtheta) by Monte Carlo
# needed to estimate covariance matrix of thetahat (regular numeric diff is not good if
# estimates are not gradient-based, and if the obj. fn. is not smooth (eg with epan kernel).
#
# usage: omega = snm_moments_jacobian(theta, data, momentargs, reps, verbose, compute_nodes)
#
# inputs:
#    theta: column vector initial parameters
#    data: data matrix
#    momentargs: (cell) additional inputs needed to compute moments.
#	      May be empty ("")
#    reps: number of replications to use (default 1000)
#    verbose: (boolean) progress monitor printed to screen (false by default)
#    compute_nodes: number of compute nodes to use in Monte Carlo draws of moments
#
# outputs:
#    jac: estimated Jacobian of average moments

function jac = snm_moments_jacobian(theta, data, momentargs, reps, verbose, compute_nodes)

	if nargin < 3 error("snm_moments_jacobian: 3 arguments required"); endif
	if nargin < 4 reps = 1000; endif # default number of reps
	if nargin < 5 verbose = false; endif
	if nargin < 6; compute_nodes = 0; endif

	n = rows(data);
	wrapper_args = {theta, n, momentargs};
	outfile = "jacobian_temp_data";
	montecarlo("snm_moments_jacobian_nodes", wrapper_args, reps, outfile, 1, verbose);
	load "jacobian_temp_data";
	shell_cmd("rm jacobian_temp_data"); # need to avoid re-use
	junk = jacobian_temp_data(:,3:columns(jacobian_temp_data));
	junk = mean(junk);
	L = momentargs{4};
	M = momentargs{6};
	jac = reshape(junk,columns(junk)/rows(theta), rows(theta));

endfunction
