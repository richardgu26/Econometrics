# ESS from linear regression
function ess = ess(y, x)
	[b, sigsq, e] = ols(y,x);
	ess = e' * e;
endfunction

# fitted values from linear regression
function fitted = fit(y,x)
	[b, sigsq, e] = ols(y,x);
	fitted = y - e;
endfunction
