## Copyright (C) 2003,2004,2005  Michael Creel <michael.creel@uab.es>
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

## usage: [V,scorecontribs,J_inv] =
##  mle_variance(theta, data, model, modelargs)
##
## This is for internal use by mle_results

# sandwich form of var-cov matrix
function [V, scorecontribs, J_inv] = mle_variance(theta, data, model, modelargs)
	scorecontribs = numgradient(model, {theta, data, modelargs});
	n = rows(scorecontribs);
	I = scorecontribs'*scorecontribs / n;
	J = numhessian("mle_obj", {theta, data, model, modelargs, 0});
	J_inv = inverse(J);
	V = J_inv*I*J_inv/n; 	% sandwich
	%V = inv(I)/n; 	% OPG
	%V = J_inv/n;     	% inv. Hessian (no minus because mle_obj is written to be minimized)
endfunction
