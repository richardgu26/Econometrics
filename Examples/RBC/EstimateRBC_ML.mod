% this is Fernández-Villaverde's basic model rbc.mod from the dynare site,
% but with psi computed to make nss = 1/3
% the other parameters are free
% I use n for labor instead if l, which looks like a 1
% also, the variable y_l (y_n) is dropped, since not used


% The file illustrates ML estimation using Dynare
% In the block "estimated_params", if you use start
% values close to the true valuse, you get reasonable
% estimates. If you try other start values, you will
% have a hard time getting convergence. This shows that
% ML may not be a good method for these models.

% true value of psi is 1.73: see ComputePsi.m

% possible explanations for difficulties
% * poor identification due to only 1 observable variable, and
%   the use of 1st order perturbation
% * ML obj function may not be nice, so hard to maximize

% Basic RBC Model 
% Michael Creel's modification of Fernández-Villaverde's original

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

% the following will lead to convergence
% don't modify this comment!
%estimated_params;
%alpha, 0.33; 
%beta, 0.99; 
%rho, 0.95; 
%sigma, 0.01;
%end;
 

% note: trying to estimate psi illustrates lack of identification
estimated_params;
alpha, 0.35; 
beta, 0.99; 
rho, 0.95;
sigma, 0.01;
%delta, 0.02;
%nss, 1/3;
%psi,  1.73;
end;

varobs c;

// computes only the posterior mode for demonstration. 
estimation(datafile=c,nobs=400, mode_compute=4, order=1) c;


