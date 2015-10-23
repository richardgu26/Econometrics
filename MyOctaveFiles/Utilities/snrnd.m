# Copyright (C) 2011 Michael Creel <michael.creel@uab.es>
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
# along with this program; If not, see <http://www.gnu.org/licenses/>.

# snrnd: random numbers from skew-normal distribution
#
# usage:
# 	z = snrnd(n, alpha, mu, sig)
#
# inputs:
#     n: number of draws to make
# alpha: SN shape parameter, positive for long R tail, negative for long L tail
#    mu: the mean of draws
#   sig: the standard error of draws
# outputs:
#     z: nx1 vector: the random draws
#
# References:
#  There is an R package for this ("sn")
#  http://azzalini.stat.unipd.it/SN/faq-r.html
#  http://azzalini.stat.unipd.it/SN/
function x = snrnd(n, alpha, m, s)
	delta = alpha / sqrt(1 + alpha^2);
	u0 = randn(n,1);
	v = randn(n,1);
	u1 = delta*u0 + sqrt(1-delta^2)*v;
	test = u0 >=0;
	z = u1.*test -u1.*(1-test);
	mu = sqrt(2/pi)*delta;
	v = (1 - 2/pi*delta^2);
	x = (z-mu)/sqrt(v);
endfunction