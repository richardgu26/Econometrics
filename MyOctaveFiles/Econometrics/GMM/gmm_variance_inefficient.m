# Copyright (C) 2003,2004  Michael Creel <michael.creel@uab.es>
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

# GMM variance, which assumes weights are not optimal
function V = gmm_variance_inefficient(theta, data, weight, omega, moments, momentargs)
	D = numgradient("average_moments", {theta, data, moments, momentargs, 0});
	D = D';
	m = feval(moments, theta, data, momentargs); # find out how many obsns. we have
	n = rows(m);
	J = D*weight*D';
	J = inv(J);
	I = D*weight*omega*weight*D';
	V = (1/n)*J*I*J;
endfunction
