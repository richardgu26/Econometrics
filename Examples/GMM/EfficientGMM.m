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
  thetahat = gmm_estimate(0, data, W, "chi2moments");
  % get the m_t at the first round estimate
  m = chi2moments(thetahat, data);
  % compute estimated covariance matrix
  % (the m_t are iid in this example)
  % and the weight matrix is the inverse
  W = inv(cov(m)); # efficient weight matrix
  % do a second round efficient estimation
  thetahat2 = gmm_estimate(0, data, W, "chi2moments");
  results = sqrt(n)*[(thetahat - theta) (thetahat2 - theta)];
endfunction  

n = 1000;
theta = 3;
args = {theta, n};
outfile = "junk.out"
reps = 1000;
system("rm junk.out"); # start from scratch
montecarlo("wrapper", args, reps, outfile, 100, false, true);
load junk.out;
hist(junk(:,3), 51);
title("Inefficient");
%print("Inefficient.png", "-dpng");
figure;
hist(junk(:,4), 51);
title("Efficient");
%print("Efficient.png", "-dpng");
