function betahat = olswrapper (args)
  n = args{1};
  theta = args{2};
  x = [ones(n,1) randn(n,1)];
  y = x*theta + randn(n,1);
  betahat = ols(y,x);
  betahat = betahat';
endfunction

