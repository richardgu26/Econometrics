# Wald test
#
# written by M. Creel, Feb. 4, 2004.  michael.creel@uab.es
#
# purpose: Wald test of restrictions Rb=r
# inputs:
# 		 b: kx1 parameter vector
# 		 varb: kxk estimated covariance matrix
# 		 R: R above, a qxk matrix
# 		 r: r above, a qx1 vector
#
# returns:

function W = WaldTest(b, varb, R, r)

	q = rows(R);

  	W = (R*b-r)'*inv(R*varb*R')*(R*b-r);

	pvalue = 1 - chi2cdf(W,q);
	tests = [W, pvalue];

	rlabel = char("Wald test");
	labels = char("Value","p-value");
	prettyprint(tests, rlabel, labels);
endfunction
