# GMM estimation for a sample from Chi^2(theta)
# compare to two method of moments estimators (see chi2mm.m)
1;

function m = chi2moments(theta, data)
   m1 = theta - data;
   m2 = theta - 0.5*(data - mean(data)).^2;
   m = [m1 m2];
endfunction   

n = 10;
theta = 3;
data = chi2rnd(theta, n, 1);

# MM 1
thetahat = mean(data);
printf("MM, v1: %f\n", thetahat);

# MM2
thetahat = 0.5*var(data);
printf("MM, v2: %f\n", thetahat);

# GMM
W = eye(2);
thetahat = gmm_estimate(0, data, W, "chi2moments");
printf("GMM: %f\n", thetahat);

