% Compile all files needed by the LBFGS wrapper.
%
% Copyright 2005-2006 Liam Stewart
% See COPYING for license.

fprintf('Compiling lbfgs mex files...\n');
fprintf('Change directory to lbfgs for this to work.\n');

mex lbfgs_mex.c list.c utils.c routines.f

fprintf('Done.\n');
