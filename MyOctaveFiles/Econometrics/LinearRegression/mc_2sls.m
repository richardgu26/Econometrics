# calculates 2SLS estimator
function [b, varb, e] = mc_2sls(y, z, zhat, names)


	xx_inv = inv(zhat'*z);
	b = xx_inv*zhat'*y;
	e = y - z*b;
	n = rows(z);
	k = columns(z);
	sigsq = e' * e / (n);
	varb = xx_inv*sigsq;
	seb = sqrt(diag(varb));
	t = b ./ seb;
	
	tss = y - mean(y);
	tss = tss' * tss;
	ess = e' * e;
	rsq = 1 - ess / tss;

	labels = char("estimate","st.err.","t-stat.","p-value");
	
	printf("\n*******************************************************\n");
	printf("2SLS estimation results\n");
	printf("Observations %d\n",n);
	printf("R-squared %f\n",rsq);
	printf("Sigma-squared %f\n",sigsq);
	printf("\n");
	p = 2 - 2*tcdf(abs(t), n - k);
	results = [b, seb, t, p];
	prettyprint(results, names, labels);

	printf("\n*******************************************************\n");
endfunction
