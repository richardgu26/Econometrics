# GMM estimation for a sample from Chi^2(theta)
# compare to two method of moments estimators (see chi2mm.m)
1;

function m = chi2moments(theta, data)
   m1 = theta - data;
   m2 = theta -  0.5*(data - mean(data)).^2;
   m = [m1 m2];
endfunction   


function results = wrapper(args)
  theta = args{1};
  n = args{2};
  data = chi2rnd(theta, n, 1);
  W = eye(2);
  results = gmm_estimate(0, data, W, "chi2moments");
  results = sqrt(n)*(results - theta);
endfunction  

n = 30;
theta = 3;
args = {theta, n};
outfile = "gmm_as_norm.out"
reps = 1000;
system("rm gmm_as_norm.out"); # start from scratch
montecarlo("wrapper", args, reps, outfile, 100, true, true);
load gmm_as_norm.out;
hist(gmm_as_norm(:,3), 51);
%if (n == 10) print ("AsNorm_n10.png", "-dpng"); endif
%if (n == 1000) print ("AsNorm_n1000.png", "-dpng"); endif
