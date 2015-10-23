# Copyright (C) 2013 Michael Creel <michael.creel@uab.es>
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

# knn_regression: k nearest neighbors regression estimator
#
# usage:
# 	[fit se min max] = knn_regression(eval_points, depvars, condvars)
#
# inputs:
#	eval_points: PxK matrix of points at which to calculate the density
#	depvars: NxG vector of observations of the G dependent variables
#	condvars: NxK matrix of data points
#	k (optional): positive scalar, number of neighbors.
#		Default is 1.5*N ^0.25  WARNING - QUITE ARBITRARY. NEEDS INVESTIGATION
#	method (optional): 1,2 or 3
#		1 = plain mean of neighbors,
#		2 = plain median of neighbors
#		3 (default) = kernel weight so closer neighbors have higher weight
#	prewhiten bool (optional): default true. If true, rotate data
# 		using Choleski decomposition of inverse of covariance,
#		to approximate independence after the transformation
#	do_cv: bool (optional). default false. If true, calculate leave-1-out
#		 fit to calculate the cross validation score
#   AIS_weights: prior/(AIS density), used to to SBIL via AIS. Default sets to 1 (not used) 
# outputs:
#	fit: PxG matrix: the fitted value at each of the P evaluation points, for \
#		each of the G dependent variables.
#	se: PxG matrix: the standard error of the neighbors at each of the P evaluation points, for \
#		each of the G dependent variables.
#	min_k: PXG matrix: the minimum of the neighbors at each of the P evaluation ponts, for \
#		each of the G dependent variable
#	max_k: PXG matrix: the maximum of the neighbors at each of the P evaluation ponts, for \
#		each of the G dependent variable
##


function [fit se] = knn_regression(eval_points, depvars, condvars, k, method=3, prewhiten=true, do_cv=false)

	if nargin < 3; error("knn_estimate: at least 3 arguments are required"); endif

	[n, G] = size(depvars);
	nn = rows(eval_points);

	# set defaults for optional args
	if (nargin < 4) k = floor(1.5*(n^(0.25))); endif	# number of neighbors
	if !isnumeric(k) k = floor(1.5*(n^(0.25))); endif # allow using "" as a placeholder

	
	if prewhiten
		[condvars scalecoefs] = scale_data([condvars; eval_points]);
		eval_points = condvars(n+1:n+nn,:);
		condvars = condvars(1:n,:);
	endif

	kk = k;
	if do_cv
		printf("doing CV, so eval_points set equal to condvars\n");
		eval_points = condvars;
		kk = k+1; % increment, because we'll drop self when doing CV
	endif


	fit = zeros(nn,G);
	se = zeros(nn,G);
	min_n = zeros(nn,G);
	max_n = zeros(nn,G);

	% the indices are in a nnXk matrix
	[nn_idx, dd] = nearest_neighbors(eval_points, condvars, kk);

	# if do_cv, we need to not use self as a neighbor
	if do_cv
		nn_idx = nn_idx(:,2:columns(nn_idx));
	endif

	% we want neighbors in nn block of k
	nn_idx = vec(nn_idx');
	dd = vec(dd');
	keepers = depvars(nn_idx,:);
	% get fits
	for i = 1:nn
		keepers_i = keepers(i*k-k+1:i*k,:);
		% plain mean
		if method == 1
			fit(i,:) = mean(keepers_i);
			se(i,:) = std(keepers_i);
		% plain median
		elseif method == 2
			fit(i,:) = median(keepers_i);
			se(i,:) = std(keepers_i);
		% Gaussian kernel to weight by distance
		elseif method == 3
			weight = dd(i*k-k+1:i*k,:);
			if max(weight) > 0
				weight = 2*weight/max(weight);
			else
			weight = ones(size(weight));
			endif
			weight = normpdf(weight);
			weight = weight/sum(weight(:));
			fit(i,:) = sum(diag(weight)*keepers_i);
			se(i,:) = std(keepers_i);
			min_k(i,:) = min(keepers_i);
			max_k(i,:) = max(keepers_i);
		endif
	endfor
endfunction
