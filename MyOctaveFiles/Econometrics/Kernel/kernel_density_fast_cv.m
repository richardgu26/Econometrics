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

# kernel_density_fast_cv: find optimal bandwith doing
#		cross validation on a randomly selected subsample
# inputs:
#	* nobs: integer - size of sample to fit to
#	* data: data matrix
#	* kernel (optional, string) the kernel function to use
# output:
#	* h: the optimal bandwidth

function bandwidth = kernel_density_fast_cv(nobs, data)

	junk = [rand(rows(data),1) data];
	junk = sortbyc(junk, 1);
	data_to_fit = junk(1:nobs,2:columns(junk));
	data = junk(nobs+1:rows(junk),2:columns(junk));
	clear junk;

	# SA controls
	ub = 3;
	lb = -5;
	nt = 1;
	ns = 1;
	rt = 0.05;
	maxevals = 50;
	neps = 5;
	functol = 1e-2;
	paramtol = 1e-3;
	verbosity = 0;
	minarg = 1;
	sa_control = { lb, ub, nt, ns, rt, maxevals, neps, functol, paramtol, verbosity, 1};

	# bfgs controls
	bfgs_control = {10};

	logbandwidth = samin("kernel_density_fast_cvobj", {1, data_to_fit, data}, sa_control);
#	logbandwidth = bfgsmin("kernel_density_fast_cvobj", {logbandwidth, data_to_fit, data}, bfgs_control);
	bandwidth = exp(logbandwidth);
endfunction
