% this is Fernández-Villaverde's basic model rbc.mod from the dynare site,
% but with psi computed to make nss = 1/3
% the other parameters are free
% I use n for labor instead if l, which looks like a 1
% also, the variable y_l (y_n) is dropped, since not used

% Basic RBC Model 
%
% Jesus Fernandez-Villaverde
% Philadelphia, March 3, 2005

%----------------------------------------------------------------
% 0. Housekeeping (close all graphic windows)
%----------------------------------------------------------------

close all;

%----------------------------------------------------------------
% 1. Defining variables
%----------------------------------------------------------------

var y c k i n z;
varexo e;

// Parameters
parameters alpha beta delta psi rho sigma nss;
alpha   = 0.33;
beta    = 0.99;
delta   = 0.023;
rho     = 0.95;  
sigma   = 0.010606;
nss = 1/3;


% this is the part that finds psi to make lss equal specified value
phi = ((1/alpha)*(1/beta -1 + delta)) ^ (1/(1-alpha));
omega = phi^(1-alpha) - delta;
kss = nss/phi;
css = omega*kss;
yss = kss^alpha * nss^(1-alpha);
psi = (1 - nss)/css * (1-alpha)*kss^alpha * nss^(-alpha);



%----------------------------------------------------------------
% 3. Model
%----------------------------------------------------------------

model; 
  (1/c) = beta*(1/c(+1))*(1+alpha*(k^(alpha-1))*(exp(z(+1))*n(+1))^(1-alpha)-delta);
  psi*c/(1-n) = (1-alpha)*(k(-1)^alpha)*(exp(z)^(1-alpha))*(n^(-alpha));
  c+i = y;
  y = (k(-1)^alpha)*(exp(z)*n)^(1-alpha);
  i = k-(1-delta)*k(-1);
  z = rho*z(-1)+sigma*e;
end;

%----------------------------------------------------------------
% 4. Computation
%----------------------------------------------------------------

initval;
  k = kss;
  c = css;
  n = nss;
  y = yss;
  z = 0; 
  e = 0;
end;

shocks;
var e = 1;
end;

steady;


estimated_params;
alpha, beta_pdf, 0.3, 0.1; 
//alpha, uniform_pdf, , ,0.3,0.5; 
beta, beta_pdf, 0.95, 0.03;
//beta, uniform_pdf, , ,0.9,0.999; 
delta, beta_pdf, 0.03, 0.01;
//delta, uniform_pdf, , ,0.02,0.06; 
rho, beta_pdf, 0.9, 0.05;
//rho, uniform_pdf, , ,0.8,0.999; 
sigma, 0.010606, inv_gamma_pdf, 0.015, inf;
nss, beta_pdf, 0.3, 0.05; 
end;

varobs c;

// short MCMC for illustration only 
estimation(datafile=c,nobs=160, order=2, mh_replic=500, mh_nblocks=1) y c k i n;


