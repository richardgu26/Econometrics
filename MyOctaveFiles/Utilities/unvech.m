## Copyright (C) 2006  Michael Creel <michael.creel@uab.es>
## Copyright (C) 2009  Jaroslav Hajek
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
## along with this program; If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} {} unvech (@var{v})
## Performs the reverse of "vech". Generates a symmetric matrix from the lower
## triangular elements, received as a vector @var{v}.
## @end deftypefn

function x = unvech (v)

  if (nargin != 1)
    usage ("unvech (v)");
  endif
  
  if (! isvector(v))
    usage ("unvech (v)");
  endif
  
  # find out dimension of symmetric matrix
  p = length (v);
  n = -(1 - sqrt (1 + 8*p))/2;
  
  if (mod (n, 1) != 0)
    error("unvech: the input vector does not generate a square matrix");
  endif
  
  x = zeros (n, n);

  # do the reverse of vech
  count = 0;
  for j = 1 : n
    i = j : n;
    x(j,i) = x(i,j) = v(count + i);
    count += n - j;
  endfor
endfunction
