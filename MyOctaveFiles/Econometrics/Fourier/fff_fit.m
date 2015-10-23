function fit = fff_fit(y, x, kalpha, evalpoints)
	if (nargin > 3)
		n = rows(x);
		nn = rows(evalpoints);
		# chop evalpoints to fit into limits of x
		m = min(x);
		m = repmat(m, nn, 1);
		test = evalpoints < m ;
		evalpoints = (1 - test) .* evalpoints + test .* m;
		m = max(x);
		m = repmat(m, nn, 1);
		test = evalpoints > m ;
		evalpoints = (1 - test) .* evalpoints + test .* m;
		x = [x; evalpoints];
	endif
	z = fff_regressors(x, kalpha);
	if (nargin > 3)
		zz = z(1:n,:);
		zfit = z(n + 1:rows(z),:);
	else
		zz = z;
		zfit = z;
	endif
	b = ols(y, zz);

	fit = zfit*b;
endfunction


