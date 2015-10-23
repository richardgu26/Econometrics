# Copyright (C) 2010 Michael Creel <michael.creel@uab.es>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA

# usage: [b, varb, e, ess] = mc_olsr(y, x, R, r, names, silent, regularvc)
# Calculates restricted LS estimator (subject to Rb=r) using the Huber-White heteroscedastic
# consistent variance estimator.
#
# inputs:
# y: dep variable
# x: matrix of regressors
# R: matrix R in Rb=r
# r: vector r in Rb=r
# names (optional) names of regressors
# silent (bool) default false. controls screen output
# regularvc (bool) default false. use normal varcov estimator, instead of het consistent (default)
#
# outputs:
# b: estimated coefficients
# varb: estimated covariance matrix of coefficients (Huber-White by default, ordinary OLS if requested with switch)
# e: ols residuals
# ess: sum of squared residuals

function [b, varb, e] = mc_olsr(y, x, R, r, names, silent, regularvc)

	if nargin < 7 regularvc = false; end
	if nargin < 6 silent = false; end
	k = columns(x);
	if (nargin < 5) || (rows(names) != k)
		names = 1:k;
		names = names';
	end

	[b, sigsq, e] = ols(y, x);

	xx_inv = inv(x'*x);
	n = rows(x);
	k = columns(x);
	q = rows(R);

	P_inv = inv(R*xx_inv*R');
	b = b - xx_inv*R'*P_inv*(R*b-r);

	e = y-x*b;
	ess = e' * e;
	sigsq = ess/(n - k - q);

	% Ordinary or het. consistent variance estimate
	if regularvc
		varb = xx_inv*sigsq;
	else
		varb = HetConsistentVariance(x,e);
	endif

	A = eye(k) - xx_inv*R'*P_inv*R;  # the matrix relating b and b_r
	varb = A*varb*A';
	seb = sqrt(diag(varb));
	t = b./seb;

	tss = y - mean(y);
	tss = tss' * tss;
	rsq = 1 - ess / tss;

	labels = char('estimate', 'st.err.', 't-stat.', 'p-value');
	if !silent
		printf('\n*********************************************************\n');
		printf('Restricted LS estimation results\n');
		printf('Observations %d\n',n);
		printf('R-squared %f\n',rsq);
		printf('Sigma-squared %f\n',sigsq);
		p = 2 - 2*tcdf(abs(t), n - k - q);
		results = [b, seb, t, p];
		if regularvc
			printf('\nResults (Ordinary var-cov estimator)\n\n');
		else
			printf('\nResults (Het. consistent var-cov estimator)\n\n');
		end
		prettyprint(results, names, labels);
		printf('\n*********************************************************\n');
	endif
endfunction
