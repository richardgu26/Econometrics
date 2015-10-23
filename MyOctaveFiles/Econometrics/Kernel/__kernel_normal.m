# Copyright (C) 2006 Michael Creel <michael.creel@uab.es>
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

# __kernel_normal: this function is for internal use by kernel_density
# and kernel_regression
#
# product normal kernel
# input: PxK matrix - P data points, each of which is in R^K
# output: Px1 vector, input matrix passed though the kernel
# other multivariate kernel functions should follow this convention

function z = __kernel_normal(z)

	z = normpdf(z);
	z = prod(z,2);

endfunction
