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


# The GMM objective function, for internal use by gmm_estimate
# This is scaled so that each piece converges to a finite number. Needs to be 
# multiplied by the sample size to get the chi-square statistic

# adapted to use lbfgs 10-May-2010. Still need to parallelize and allow
# moments to return analytic derivative

function [obj_value, score] = gmm_obj(theta, data, weight, moments, momentargs, compute_nodes)
	m = average_moments(theta, data, moments, momentargs, compute_nodes);
	obj_value = m' * weight *m;
	score = numgradient("average_moments", {theta, data, moments, momentargs, compute_nodes});
	score = score'*weight*m;

	if (((abs(obj_value) == Inf)) || (isnan(obj_value)))
		obj_value = realmax;
	endif	
endfunction	
