function [stop,user_data] = lbfgs_null_cb(x,iter,state,user_data,opts)
% LBFGS_NULL_CB  Callback that does nothing.
%
% [stop,ud] = lbfgs_null_cb(x,iter,state,ud,opts)
%
% Copyright 2005-2006 Liam Stewart
% See COPYING for license.

stop = 0;
user_data = [];
