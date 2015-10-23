## Copyright (C) 2013 Michael Creel
## 
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with Octave; see the file COPYING.  If not, see
## <http://www.gnu.org/licenses/>.

## ReadBinaryData
##
## usage: data = ReadBinaryData('filename');
##
## reads data written by WriteBinaryData. The first number
## in the data set is the record size, p, and then there are
## some number, n, of records of that size. This returns the
## data in a nXp matrix.

## Author: Michael Creel <michael@pelican>
## Created: 2013-09-04

function [ data ] = ReadBinaryData(filename)
	FN=fopen(filename, 'r');
	ncols = fread(FN, 1, 'double');
	data = fread(FN, Inf);
	fclose(FN);
	data = reshape(data, rows(data)/ncols, ncols);
	data = data;
endfunction
