## Copyright (C) 2003,2004,2005  Michael Creel <michael.creel@uab.es>
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
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

# usage:
# [theta, obj_value, conv, iters] = mle_estimate(theta, data, model, modelargs, control, nslaves)
#
# inputs:
# theta: column vector of model parameters
# data: data matrix
# model: name of function that computes log-likelihood
# modelargs: (cell) additional inputs needed by model. May be empty ("")
# control: (optional) BFGS or SA controls (see bfgsmin and samin). May be empty ("").
# nslaves: (optional) number of slaves if executed in parallel (requires MPITB)
#
# outputs:
# theta: ML estimated value of parameters
# obj_value: the value of the log likelihood function at ML estimate
# conv: return code from bfgsmin (1 means success, see bfgsmin for details)
# iters: number of BFGS iteration used
#
# please see mle_example.m for examples of how to use this
function [theta, obj_value, convergence, iters] = mle_estimate(theta, data, model, modelargs, control, nslaves)


	if nargin < 3
		error("mle_estimate: 3 arguments required");
	endif

	if nargin < 4 modelargs = {}; endif # create placeholder if not used
	if !iscell(modelargs) modelargs = {}; endif # default controls if receive placeholder
	if nargin < 5 control = {-1,0,1,1}; endif # default controls and method
	if !iscell(control) control = {-1,0,1,1}; endif # default controls if receive placeholder
	if nargin < 6 nslaves = 0; endif
	if nslaves > 0
		global NSLAVES PARALLEL NEWORLD TAG;
		LAM_Init(nslaves);
		# Send the data to all nodes
		NumCmds_Send({"data", "model", "modelargs"}, {data, model, modelargs});
	endif

	# bfgs or sa?
	if (size(control,1)*size(control,2) == 0) # use default bfgs if no control
		control = {Inf,0,1,1};
		method = "bfgs";
	elseif (size(control,1)*size(control,2) < 11)
		method = "bfgs";
	else method = "sa";
	endif

	# do estimation using either bfgsmin or samin
 	if strcmp(method, "bfgs")
		% to use bfgsmin
		%[theta, obj_value, convergence, iters] = bfgsmin("mle_obj", {theta, data, model, modelargs, nslaves}, control);
		
    % to use fminunc
    options = optimset('TolFun', 1e-8);
		options = optimset(options, 'TolX', 1e-5);
		options = optimset(options, 'AutoScaling', 'on');
		[theta, obj_value, convergence, output] = fminunc(@(theta) mle_obj(theta, data, model, modelargs, nslaves), theta, options);
	  convergence = convergence > 0;
    
    % to use lbfgs
    % maxits = control{1};		
	 	%opt_fun = @(t) mle_obj (t, data, model, modelargs, nslaves);
	  %opts = lbfgs_options('iprint', -1, 'maxits', maxits, 'factr', 1e10, 'cb', @mle_callback);
 	  %	[theta, obj_value, convergence, userdata] = lbfgs(opt_fun, theta, opts);	
		%if (convergence == 0) convergence = 1; endif
		%ters = userdata.its;
	elseif strcmp(method, "sa")
		[theta, obj_value, convergence] = samin("mle_obj", {theta, data, model, modelargs, nslaves}, control);
	endif
	


	if nslaves > 0
		LAM_Finalize;
	endif # cleanup
	obj_value = - obj_value; # recover from minimization rather than maximization
endfunction
