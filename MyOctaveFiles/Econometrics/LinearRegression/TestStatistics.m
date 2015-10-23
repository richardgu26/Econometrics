# TESTS - the classical test statistics for linear restrictions
# 
# written by M. Creel, March 10, 1998.  michael.creel@uab.es
# converted to Octave 13 Jan. 2003
# 
# purpose: calculates F, Wald, Score and Likelihood Ratio tests for
# linear model						 y=XB+e
# 									 e~N(0,sig^2*I_n)
# subject to linear restrictions  	 RB=r
# 
# format: [F, W, LR, S] = TestStatistics(y,x,R,r);
# 
# inputs:
# 		 y: nx1 dependent variable
# 		 x: nxk regressor matrix
# 		 R: R above, a qxk matrix
# 		 r: r above, a qx1 vector
# 
# returns: F: the F statistic
# 		 W: the Wald statistic
# 		 S: the score statistic
# 		 LR: the likelihood ratio statistic

function [F, W, LR, S] = TestStatistics(y, x, R, r)

    n = rows(x);
    k = columns(x);
    q = rows(R);

    # OLS
	[b, sigsq, e] = ols(y, x);

	# The restricted estimator
	xx_inv = inv(x'*x);
	P_inv = inv(R*xx_inv*R');
	b_r = b - xx_inv*R'*P_inv*(R*b-r);

  # Sums of squared errors and estimators of sig^2
  e = y - x*b;
  ess = e'*e;
  e_r = y - x*b_r;
  ess_r = e_r' * e_r;
  sigsqhat_ols = ess/(n-k);
  sigsqhat_mle = ess/(n);
  sigsqhat_mle_r = ess_r/(n);

  # F-test
  F = (ess_r-ess)/q;
  F = F/sigsqhat_ols;

  # Wald test (uses unrestricted model's est. of sig^2
  W = (R*b-r)'*P_inv*(R*b-r)/sigsqhat_mle;

  # Score test (uses restricted model's est. of sig^2 
  P_x = x * xx_inv * x';
  S = e_r' * P_x * e_r/(sigsqhat_mle_r);

  # LR test
  lnl = -n/2*log(2*pi) - n/2*log(sigsqhat_mle) - e' * e/(2*sigsqhat_mle);
  lnl_r = -n/2*log(2*pi) - n/2*log(sigsqhat_mle_r) - e_r' * e_r/(2*sigsqhat_mle_r);
  LR = 2*(lnl-lnl_r);
	tests = [F;W;LR;S];
	WLRS = [W;LR;S];
	pvalues = [1 - fcdf(F,q,n-k); 1 - chi2cdf(WLRS,q)]; 
	tests = [tests, pvalues];
	
	TESTS = char("F","Wald","LR","Score");
	labels = char("Value","p-value");
	prettyprint(tests, TESTS, labels);
endfunction
