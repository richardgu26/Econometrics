function [x,f,exitflag,user_data] = lbfgs(fn, x0, lbd, ubd, nbd, opts)
% LBFGS  Minimize a function using the L-BFGS algorithm.
%
% [x,fx,ef,ud] = LBFGS(fn, x0, opts)
% [x,fx,ef,ud] = LBFGS(fn, x0, lb, ub, nbd, opts)
%
% fn   - handle to function to minimize
% x0   - starting point
% lb   - vector of lower bounds
% ub   - vector of upper bounds
% nbd  - vector of bound types
%          nbd(i) = 0 if x(i) is unbounded
%                   1 if x(i) has a lower bound (in lb(i))
%                   2 if x(i) has both lower (lb(i)) and upper (ub(i)) bounds
%                   3 if x(i) has an upper bound (ub(i))
% opts - L-BFGS option structure (see LBFGS_OPTIONS) (optional)
%
% x    - final value
% fx   - fn(x)
% ef   - exit flag.
%          0: algorithm converged or max iterations reached
%          1: abnormal termination
%          2: error
%          3: unknown cause of termination
% ud   - user data from callback function
%
% Copyright 2005-2006 Liam Stewart
% See COPYING for license.

if nargin < 6; opts = lbfgs_options; end;
if nargin < 5; nbd = []; end;
if nargin < 4; ubd = []; end;
if nargin < 3; lbd = []; end;
if nargin < 2; error('not enough arguments'); end;

if isstruct(lbd) % backwards compatility
    opts = lbd;
    lbd = [];
end

n = length(x0(:));
f = 0;
g = zeros(n, 1);

if isempty(nbd)
    nbd = zeros(n,1);
end
if isempty(ubd)
    ubd = zeros(n,1);
end
if isempty(lbd)
    lbd = zeros(n,1);
end

assert(length(nbd) == n);
assert(length(ubd) == n);
assert(length(lbd) == n);

id = lbfgs_mex('init', n, x0, lbd, ubd, nbd, opts);

iter = struct('it', 0, 'f', f, 'g', g);
[stop,user_data] = opts.cb(x0,iter,'init',[],opts);

maxits = opts.maxits;
exitflag = 0;
it = 1;
while 1
    [x,cmd] = lbfgs_mex('step', id, f, g);
    
    if cmd == 0                         % evaluate f, g
        [f,g] = fn(x);
    elseif cmd == 1                     % iteration complete
        iter = struct('it', it, 'f', f, 'g', g);
        [stop,user_data] = opts.cb(x,iter,'iter',user_data,opts);

        if it >= maxits || stop
            exitflag = 0;
            lbfgs_mex('stop', id);
            break;
        else
            it = it + 1;
        end
    elseif cmd == 2                     % converged
        exitflag = 0;
        break;
    elseif cmd == 3                     % abnormal termination
        exitflag = 1;
        break;
    elseif cmd == 4                     % error
        exitflag = 2;
        break;
    else                                % unknown
        exitflag = 3;
        break;
    end
end

iter = struct('it', it, 'f', f, 'g', g);
[stop,user_data] = opts.cb(x,iter,'done',user_data,opts);

lbfgs_mex('destroy', id);
