// Copyright (C) 2007  Michael Creel <michael.creel@uab.es>
// under the terms of the GNU General Public License.
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; If not, see <http://www.gnu.org/licenses/>.

// __kernel_weights: for internal use by kernel density and regression functions

#include <oct.h>
#include <octave/parse.h>

DEFUN_DLD(__kernel_weights, args, ,"__kernel_weights: for internal use by kernel_regression and kernel_density functions")
{

	int nargin = args.length();
	Matrix data (args(0).matrix_value());
	Matrix evalpoints (args(1).matrix_value());
	std::string kernel (args(2).string_value());
        int startblock = 0;
        // startblock > 0 indicates that leave-1-out fit should be used
        if (nargin > 3) {
            int startblock (args(3).int_value());
        }    

	int n, nn, i, j, k, kk, dim;

	n = data.rows();
	dim = data.columns();
	nn = evalpoints.rows();
	Matrix W(nn, n);
	Matrix zz(n,dim);
	Matrix zzz;
	octave_value kernelargs;
	octave_value_list f_return;

	for (i = 0; i < nn; i++) {
		for (j = 0; j < n; j++) {
			// note to self: zz.insert(data.row(j) - evalpoints.row(i), j, 0) is slower
			for (k = 0; k < dim; k++) {
				zz(j,k) = data(j,k) - evalpoints(i,k);
			}
		}
		kernelargs = zz;
		f_return = feval(kernel, kernelargs);
		zzz = f_return(0).matrix_value();
		// note to self: W.insert(zzz.transpose(),i,0) is slower
		for (j = 0; j < n; j++) {
			W(i,j) = zzz(j,0);
		}
                if (startblock > 0) W(i,i+startblock-1) = 0.0; // leave self out of fit
	}
	f_return(0) = W;
	return f_return;
}

