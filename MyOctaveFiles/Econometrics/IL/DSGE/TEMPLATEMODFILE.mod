// The two shock model of the paper
// This takes steady state of hours as given, and computes the psi
// parameter that corresponds.
close all;

// Define variables
var y c k n i z1 z2 MUC MUL MPK MPL;
varexo e1 e2;

// Parameters
parameters alpha beta delta gam  nss rho1 sigma1 rho2 sigma2 psi c1 iss yss kss css;
load parameterfile1;
set_param_value('alpha',alpha)
set_param_value('beta',beta)
set_param_value('delta',delta)
set_param_value('gam',gam)
set_param_value('rho1',rho1)
set_param_value('sigma1',sigma1)
set_param_value('rho2',rho2)
set_param_value('sigma2',sigma2)
set_param_value('nss',nss)

// find the psi that makes steady state labor (nss) equal to the specified value
c1 = ((1/beta + delta - 1)/alpha)^(1/(1-alpha));
kss = nss/c1;
iss = delta*kss;
yss = kss^alpha * nss^(1-alpha);
css = yss - iss;
psi =  (css^(-gam)) * (1-alpha) * (kss^alpha) * (nss^(-alpha));

// Model
model;
  MUC = (c)^(-gam);
  MUL = psi*exp(z2);
  MPK = alpha * exp(z1) * (k(-1))^(alpha-1) * n^(1-alpha);  
  MPL = (1-alpha)*exp(z1)* (k(-1))^alpha * n^(-alpha); 
  MUC = beta*MUC(+1) * (1 + MPK(+1) - delta);
  MUL/MUC = MPL;
  z1 = rho1*z1(-1) + sigma1*e1;
  z2 = rho2*z2(-1) + sigma2*e2;
  y = exp(z1) * ((k(-1))^alpha) * (n^(1-alpha));
  i = y - c;
  k = i + (1-delta)*k(-1);
end;

// Initial values
initval;
k = kss;
n = nss;
c = css;
y = yss;
i = iss;
MUC = c^(-gam);
MPK = alpha  * k^(alpha-1) * n^(1-alpha);  
MPL = (1-alpha)* (k)^alpha * n^(-alpha); 
MUL = MPL*MUC;
z1 = 0;
z2 = 0;
end;

// Shocks
shocks;
var e1 = 1;
var e2 = 1;
end;

steady(maxit=1000, solve_algo=0);

// keep rng sequences in order on each node
load state1;
ss = round(s(1)/s(2)*sum(s));
set_dynare_seed(ss);

// generate data
stoch_simul(nograph, noprint, order=3, drop=100, periods=260) y c n MPK MPL;
