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
 psi

