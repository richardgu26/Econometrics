# Copyright (C) 2006, 2010  Michael Creel <michael.creel@uab.es>
# under the terms of the GNU General Public License.
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

# kernel_density_nodes: for internal use by kernel_density - does calculations on nodes

function z = kernel_density_nodes(myeval, data, kernel, startblock)
	# startblock coming in as arg is signal to do leave-1-out)
	if (nargin < 4)
  		W = __kernel_weights(data, myeval, kernel);
		do_cv = false;
	else
  		W = __kernel_weights(data, myeval, kernel, startblock);
		do_cv = true;
	endif		  
	z = sum(W,2);
	n = rows(data);
	if do_cv
		z = z/(n-1);
	else
		z = z/n;
	endif
endfunction


