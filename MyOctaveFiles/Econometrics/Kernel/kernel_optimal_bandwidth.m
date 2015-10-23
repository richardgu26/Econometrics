# Copyright (C) 2007 Michael Creel <michael.creel@uab.es>
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

# kernel_optimal_bandwidth: find optimal bandwith doing leave-one-out cross validation
# inputs:
#	* data: data matrix
#	* depvar: column vector or empty ("").
#		If empty, do kernel density, orherwise, kernel regression
#	* kernel (optional, string) the kernel function to use
# output:
#	* h: the optimal bandwidth

function bandwidth = kernel_optimal_bandwidth(data, depvar, kernel, prewhiten)

	if (nargin < 2) error("kernel_optimal_bandwidth: 2 arguments required"); endif
	if (nargin < 3) kernel = "__kernel_normal"; endif
	if (nargin < 4) prewhiten = false; endif
	
	do_density = false;
	if isempty(depvar) do_density = true; endif;

	# SA controls
	ub = 3;
	lb = -10;
	nt = 1;
	ns = 1;
	rt = 0.05;
	maxevals = 50;
	neps = 5;
	functol = 1e-2;
	paramtol = 1e-3;
	verbosity = 2;
	minarg = 1;
	sa_control = { lb, ub, nt, ns, rt, maxevals, neps, functol, paramtol, verbosity, 1};

	# bfgs controls
	bfgs_control = {10};

	if do_density
		bandwidth = samin("kernel_density_cvscore", {1, data, kernel, prewhiten}, sa_control);
	else
		bandwidth = samin("kernel_regression_cvscore", {1, data, depvar, kernel, prewhiten}, sa_control);
	endif
	bandwidth = exp(bandwidth);
endfunction
