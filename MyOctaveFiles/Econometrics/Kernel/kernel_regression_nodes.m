# Copyright (C) 2006, 2007, 2009  Michael Creel <michael.creel@uab.es>
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

# kernel_regression_nodes: for internal use by kernel_regression - does calculations on nodes

function fit = kernel_regression_nodes(myeval, depvar, condvars, kernel, startblock)
	# startblock coming in as arg is signal to do leave-1-ou)
	if (nargin < 5)
  		W = __kernel_weights(condvars, myeval, kernel);
	else
  		W = __kernel_weights(condvars, myeval, kernel, startblock);
	endif		  
	den = sum(W,2);
	if !all(den)
		warning("kernel_regression: some evaluation points have no neighbors - increase the bandwidth");
		den = den + eps; # avoid divide by zero
	endif
	W = diag(1 ./ den)*W;
	fit = W*depvar;
endfunction


