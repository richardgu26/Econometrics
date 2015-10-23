function opts = lbfgs_options(varargin)
% LBFGS_OPTIONS  Options for controlling LBFGS.
%
% opts = LBFGS_OPTIONS('param1', val1, 'param2', val2, ...)
%
% Parameters include:
%   maxits - maximum number of iterations (default: 100)
%   m      - number of corrections to use (default: 5)
%   factr  - tol. of function value termination test (default: 1e7).
%   pgtol  - tol. of gradient projection termination test (default: 1e-5)
%   iprint - control on frequency and type of output (default: 0)
%   cb     - callback function (default: lbfgs_null_cb)
%
% See L-BFGS documentation for details on m, factr, pgtol, and iprint.
%
% Copyright 2005-2006 Liam Stewart
% See COPYING for license.

[maxits,m,factr,pgtol,iprint,cb] = ...
    process_options(varargin, ...
                    'maxits', 100, ...
                    'm', 5, ...
                    'factr', 1e7, ...
                    'pgtol', 1e-5, ...
                    'iprint', 0, ...
                    'cb', @lbfgs_null_cb);

opts.maxits = maxits;
opts.m = m;
opts.factr = factr;
opts.pgtol = pgtol;
opts.iprint = iprint;
opts.cb = cb;
