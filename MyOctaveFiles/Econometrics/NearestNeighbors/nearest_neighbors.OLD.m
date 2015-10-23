# Copyright (C) 2012 Michael Creel <michael.creel@uab.es>
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

# nearest_neighbors: find K nearest neighbors to haystack in needles
#
# usage:
# 	[idx, dist] = nearest_neighbors(needle, haystack, neighbors)
#
# inputs:
#	needles: mxK matrix: m points in R^K
#	haystack: nxK matrix: n points in R^K
#	neighbors: integer, the number of neighbors to find
# outputs:
#	idx: m x neighbors matrix; indices of the closest neighbors points in haystack
#		to each point in needles
#   	dist: m x neighbors matrix: the distances of the closest neighbors points in haystack
#		to each point in needles

function [idx, dist] = nearest_neighbors(needles, haystack, neighbors)
	dim = columns(needles);
	limit = rows(haystack);
	# save needles and haystack to file
	save -ascii /tmp/_nearest_neighbors_needles needles;
	save -ascii /tmp/_nearest_neighbors_haystack haystack;
	system("sync");
	s = sprintf("ann_sample -d %d -max %d -nn %d -df /tmp/_nearest_neighbors_haystack -qf /tmp/_nearest_neighbors_needles", dim, limit, neighbors);
	[junk output] = system(s);	# the call to ANN
	# manipulate to get idx
	output = str2num(output);
	idx = output(:,1);
	# recover from base-0 indexing
	idx = idx + 1;
	idx = reshape(idx, neighbors, rows(idx)/neighbors);
	idx = idx';
	dist = output(:,2);
	dist = reshape(dist, neighbors, rows(dist)/neighbors);
	dist = dist';
endfunction
