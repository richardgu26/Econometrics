1;
# illustrates effects of collinearity
# try setting rho=0.1 or rho=0.9, and observe
# the effect on standard dev of coefficients
function results = wrapper(args)
  rho = args{1};
  n = args{2};
  V = [1 rho; rho 1];
  P = chol(V);
  x = randn(n,2)*P;
  # cor(x) # view correlation matrix
  x = [ones(n,1) x];
  b = ones(3,1);
  y = x*b + randn(n,1);
  b1 = mc_ols(y,x, "", true);

  z = x(:,1:2);
  b2 = mc_ols(y,z, "", true);

  R = [0 1 -1];
  r = 0;
  b3 = mc_olsr(y, x, R, r, "", true);
  results = [b1' b2' b3'] ;
endfunction  

reps = 1000;
outfile = "collinearity_mc.out";
rho = 0.9;
n = 30;
args = {rho, n};
system("rm collinearity_mc.out");
montecarlo("wrapper", args, reps, outfile, 500, true);
load collinearity_mc.out;
data = collinearity_mc(:,3:10);
printf("\ncorrelation between x2 and x3: %f\n\n",rho);
printf("descriptive statistics for 1000 OLS replications\n");
dstats(data(:,1:3));
printf("descriptive statistics for 1000 OLS replications, dropping x3\n");
dstats(data(:,4:5));
printf("descriptive statistics for 1000 Restricted OLS replications, b2=b3\n");
dstats(data(:,6:8));
hist(data(:,2),30);
print("collin_ols.png","-dpng");
figure;
hist(data(:,5),30);
print("collin_drop.png","-dpng");
figure;
hist(data(:,7),30);
print("collin_rls.png","-dpng");

