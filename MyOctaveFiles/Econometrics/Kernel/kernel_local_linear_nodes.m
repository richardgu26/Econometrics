# Copyright (C) 2010  Michael Creel <michael.creel@uab.es>
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

# kernel_local_linear_nodes: for internal use by kernel_regression - does calculations on nodes

function fit = kernel_local_linear_nodes(myeval, depvars, condvars, do_cv, kernel)
	W = __kernel_weights(condvars, myeval, kernel);
	# drop own weight for CV
	if (do_cv) W = W - diag(diag(W)); endif
	den = sum(W,2);
	if !all(den)
		warning("kernel_local_linear: some evaluation points have no neighbors - increase the bandwidth");
		den = den + eps; # avoid divide by zero
	endif
	W = diag(1 ./ den)*W;
	n = rows(condvars);
	nn = rows(myeval);
	G = columns(depvars);
	fit = zeros(nn,G);
	X = [ones(n,1) condvars];

	# loop over eval points, doing GLS with kernel weights
	for i = 1:nn
		XX = diag(W(i,:))*X;
		b = inv(X'*XX)*XX'*depvars;
		f = [1 myeval(i,:)]*b;
		fit(i,:) = f;
	endfor
endfunction


