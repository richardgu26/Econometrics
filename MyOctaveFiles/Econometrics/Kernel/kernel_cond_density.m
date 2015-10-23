# Copyright (C) 2007 Michael Creel <michael.creel@uab.es>
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

# kernel_cond_density: multivariate kernel conditional density estimator
#
# usage:
# 	dens = kernel_cond_density(eval_points, data, condvars, bandwidth)
#
#
# References:
# Wand, M.P. and Jones, M.C. (1995), 'Kernel smoothing'.
# http://www.xplore-stat.de/ebooks/scripts/spm/html/spmhtmlframe73.html

function z = kernel_cond_density(eval_points, data, condvars, bandwidth_num, bandwidth_den, kernel, prewhiten, do_cv, nodes)

	if nargin < 6 kernel = "__kernel_normal"; endif
	if (nargin < 7) prewhiten = false; endif 	# automatic prewhitening?
	if (nargin < 8)	do_cv = false; endif 		# ordinary or leave-1-out
	if (nargin < 9)	computenodes = 0; endif		# parallel?


	# joint density of all vbls
	j_density = kernel_density(eval_points, data, bandwidth_num, kernel, false, false, computenodes);

	# marginal density of conditioning variables
	m_data = data(:,condvars);
	m_eval = eval_points(:,condvars);
	m_density = kernel_density(m_eval, m_data, bandwidth_den, kernel, false, false, computenodes);

	# conditional density
	check = (m_density == 0); 	# avoid div. by zero
	m_density = m_density + check;
	z = j_density ./ m_density; 	# cond = joint / marg
	z = z .* (1 - check); 		# cond density is zero if marginal is zero
endfunction
