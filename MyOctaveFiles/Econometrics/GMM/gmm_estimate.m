# Copyright (C) 2003,2004,2005 Michael Creel <michael.creel@uab.es>
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

# usage: [theta, obj_value, convergence, iters] =
#           gmm_estimate(theta, data, weight, moments, momentargs, control, compute_nodes)
#
# inputs:
#      theta: column vector initial parameters
#       data: data matrix
#     weight: the GMM weight matrix
#    moments: name of function computes the moments
#	      (should return nXg matrix of contributions)
# momentargs: (cell, optional) additional inputs needed to compute moments.
# 	      May be empty ("")
#    control: (optional) BFGS or SA controls (see bfgsmin and samin).
#             May be empty ("").
#    compute_nodes: (optional) number of slaves if executed in parallel
#             (requires MPITB)
#
# outputs:
# theta: GMM estimate of parameters
# obj_value: the value of the gmm obj. function
# convergence: return code from bfgsmin
#              (1 means success, see bfgsmin for details)
# iters: number of BFGS iteration used
#
# please type "gmm_example" while in octave to see an example

# call the minimizing routine
function [theta, obj_value, convergence, iters] = gmm_estimate(theta, data, weight, moments, momentargs="", control={-1}, compute_nodes=0)

	if nargin < 4 error("gmm_estimate: 4 arguments required"); endif
	
 	if compute_nodes > 0
		global NSLAVES PARALLEL NEWORLD NSLAVES TAG;
		LAM_Init(compute_nodes);
		# Send the data to all nodes
		NumCmds_Send({"data", "weight", "moments", "momentargs"}, {data, weight, moments, momentargs});
	endif

	# bfgs or sa?
	if (size(control,1)*size(control,2) == 0) # use default bfgs if no control
		control = {Inf,0,1,1};
		method = "bfgs";
	elseif (size(control,1)*size(control,2) < 11)
		method = "bfgs";
	else method = "sa";
	endif

	if strcmp(method, "bfgs")
	  % to use bfgsmin
		%[theta, obj_value, convergence, iters] = bfgsmin("gmm_obj", {theta, data, weight, moments, momentargs, compute_nodes}, control);
		% to use fminunc
		%options = optimset('TolFun', 1e-8);
		%options = optimset(options, 'TolX', 1e-5);
		%options = optimset(options, 'AutoScaling', 'on');
		%[theta, obj_value, convergence, output] = fminunc(@(theta) gmm_obj(theta, data, weight, moments, momentargs, compute_nodes), theta, options);
		% to use lbfgs 	  
		opt_fun = @(t) gmm_obj (t, data, weight, moments, momentargs, compute_nodes);
		opts = lbfgs_options('iprint', -1, 'maxits', 1000, 'cb', @gmm_callback);
		[theta, obj_value, convergence, userdata] = lbfgs(opt_fun, theta, opts);	
		iters = userdata.its;
		if (convergence == 0) convergence = 1; endif
	elseif strcmp(method, "sa")
	  	[theta, obj_value, convergence] = samin("gmm_obj", {theta, data, weight, moments, momentargs, compute_nodes}, control);
	endif

	if compute_nodes > 0 LAM_Finalize; endif # clean up

endfunction
