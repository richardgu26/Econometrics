# Copyright (C) 2012,2013 Michael Creel <michael.creel@uab.es>
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
# last revision: the temporary files have random names, so that this is safe to use in
# parallel on multiple nodes.
function [idx, dist] = nearest_neighbors(needles, haystack, neighbors)
	dim = columns(needles);
	limit = rows(haystack);
	# save needles and haystack to file
	symbols = ['a':'z' 'A':'Z' '0':'9'];
	MAX_ST_LENGTH = 50;
	stLength = randi(MAX_ST_LENGTH);
	targetfile = randi(numel(symbols),[1 MAX_ST_LENGTH]);
	targetfile = symbols (targetfile);
	stLength = randi(MAX_ST_LENGTH);
	queryfile = randi(numel(symbols),[1 MAX_ST_LENGTH]);
	queryfile = symbols (queryfile);
	queryfile = strcat('/tmp/',queryfile);
	targetfile = strcat('/tmp/',targetfile);
	save('-ascii', queryfile, 'needles');
	save('-ascii', targetfile, 'haystack');
	%system("sync");
	s = sprintf("ann_sample -d %d -max %d -nn %d -df %s -qf %s", dim, limit, neighbors, targetfile, queryfile);
	[junk output] = system(s);	# the call to ANN
	
	% get rid of the files
	s = sprintf("rm %s", targetfile);
	system(s);
	s = sprintf("rm %s", queryfile);
	system(s);
	%system("sync");
	
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
