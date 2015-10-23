# R^2 from linear regression
function rsq = rsquare(y, x)
	[b, sigsq, e] = ols(y,x);

	tss = y - mean(y);
	tss = tss' * tss;
	ess = e' * e;
	rsq = 1 - ess / tss;
endfunction
