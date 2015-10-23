function Z = aux_stat(data)
  # variables coming in are  y l c k w r
  n = rows(data);
  y = data(1:n-1,1);
  l = data(1:n-1,2);
  c = data(1:n-1,3);
  k = data(1:n-1,4);
  w = data(1:n-1,5);
  r = data(1:n-1,6);
  clead = data(2:n,3);
  ylead = data(2:n,1);
  klead = data(2:n,4);
  rlead = data(2:n,6);
  
  # delta
  delta = mean((w.*l + r.*k + k - c -klead)./k);

  # beta
  beta = mean(clead ./ (c + c.*rlead - c*delta));
 
  # psi
  psi = mean((1-l).*w./c);

  # alpha rho and sig
  depvar = slog(y)-slog(l);
  regs = [slog(k)-slog(l)];
  [alpha, V, e] = ols(depvar, regs);
  rho = ols(e(2:rows(e),:),e(1:rows(e)-1,:));
  rho2 = corr(c, clead);
  rho3 = corr(y, ylead); # this is new since 2011 version of paper
#  sig = sqrt(V);

  # cost function
  c = w.*l+r.*k;
  regs = [y w r];
  mc = ols(c, regs);
  mc = mc(1,:);
  eps = (1-mc); # normally woul use inverse, but don't want outliers

  # alpha and eps
  #eps1 = mean(r.*k./y/alpha);
  #eps2 = mean(w.*l./y/(1-alpha));
  sig = std(data)';

  # assemble the aux stat
  Z = [alpha; delta; beta; psi; rho; rho2; rho3; sig; eps]; 
endfunction
