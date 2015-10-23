# Copyright (C) 2005  Michael Creel michael.creel@uab.es
# under the terms of the GNU General Public License.
# The GPL license is in the file COPYING
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
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
# this prints matrices with row and column labels

# sortbyc: sort a matrix by a given column, keeping row elements together 
# usage: b = sortbyc(a,col,mode)
# example:
# a = [1 3; 2 2; 3 1];
# [ a sortbyc(a,2)]
# ans =
# 
#   1  3  3  1
#   2  2  2  2
#   3  1  1  3



function a = sortbyc(a,col,mode);
	if nargin < 3; mode = "ascend"; endif
	[junk, i] = sort(a, mode);
	a = a(i(:,col),:);
endfunction
	
 
