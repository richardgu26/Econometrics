# Calculates generalized IV estimator using the Huber-White heteroscedastic
# consistent variance estimator.

function [b, varb, e, ess] = mc_giv(y, x, w, names, silent, regularvc)
	k = columns(x);

	if nargin < 6 regularvc = 0; endif
	if nargin < 5 silent = 0; endif
	if (nargin < 4) || (rows(names) != k)
		names = 1:k;
		names = names';
	endif

	P = w*pinv(w'*w)*w';
	yy = P*y;
	xx = P*x;

	b = ols(yy,xx);
	e = y - x*b;

	xx_inv = inv(x'*x);
	n = rows(x);

	ess = e' * e;
	sigsq = ess/n;

	# Ordinary or het. consistent variance estimate
	xw = x'*w;
	ww_inv = pinv(w'*w);

	if regularvc
		varb = inverse(xw*ww_inv*xw')*sigsq;
	else

		E = e .^2;
		E = diag(E);
		A = inverse(xw*ww_inv*xw');
		varb = A*xw*ww_inv*w'*E*w*ww_inv*xw'*A;
	endif

	seb = sqrt(diag(varb));
	t = b ./ seb;

	tss = y - mean(y);
	tss = tss' * tss;
	rsq = 1 - ess / tss;

	labels = char("estimate","st.err.", "t-stat.", "p-value");
	if !silent
		printf("\n*********************************************************\n");
		printf("Generalized IV estimation results\n");
		printf("Observations %d\n",n);
		printf("R-squared %f\n",rsq);
		printf("Sigma-squared %f\n",sigsq);
		p = 2 - 2*tcdf(abs(t), n - k);
		results = [b, seb, t, p];
		if regularvc
			printf("\nResults (Ordinary var-cov estimator)\n\n");
		else
			printf("\nResults (Het. consistent var-cov estimator)\n\n");
		endif
		prettyprint(results, names, labels);
		printf("\n*********************************************************\n");
	endif
endfunction
