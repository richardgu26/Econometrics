# Copyright (C) 2009  Michael Creel <michael.creel@uab.es>
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


# SNM estimated variance of parameters, weight matrix may be optimal or not
function V = snm_variance(theta, data, weight, moments, momentargs, reps, compute_nodes, momentcov)
	D = snm_moments_jacobian(theta, data, momentargs, reps, false, compute_nodes);

	# drop wws, keep real parameters
	K = momentargs{4};
	k = rows(theta);
	k = 1:(k-K);
	D = D(:,k);

	D = D';
	n = rows(data);
	if (nargin < 8)
		V = (1/n)*inv(D*weight*D');
	else
		J = D*weight*D';
		J = inv(J);
		I = D*weight*momentcov*weight*D';
		V = (1/n)*J*I*J;
	endif
endfunction
