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

# usage: [b, varb, e, ess] = mc_ols(y, x, names, silent, regularvc)
# Calculates ordinary LS estimator using the Huber-White heteroscedastic
# consistent variance estimator.
#
# inputs:
# y: dep variable
# x: matrix of regressors
# names (optional) names of regressors
# silent (bool) default false. controls screen output
# regularvc (bool) default false. use normal varcov estimator, instead of het consistent (default)
#
# outputs:
# b: estimated coefficients
# varb: estimated covariance matrix of coefficients (Huber-White by default, ordinary OLS if requested with switch)
# e: ols residuals
# ess: sum of squared residuals

function [b, varb, e, ess, rsq] = mc_ols(y, x, names, silent, regularvc)
	k = columns(x);

	if nargin < 5 regularvc = 0; end
	if nargin < 4 silent = 0; end
	if (nargin < 3) || (rows(names) != k)
		names = 1:k;
		names = names';
	end

	[b, sigsq, e] = ols(y,x);

	xx_inv = inv(x'*x);
	n = rows(x);

	ess = e' * e;

	% Ordinary or het. consistent variance estimate
	if regularvc
		varb = xx_inv*sigsq;
	else
		varb = HetConsistentVariance(x,e);
	endif

	seb = sqrt(diag(varb));
	t = b ./ seb;

	tss = y - mean(y);
	tss = tss' * tss;
	rsq = 1 - ess / tss;

	labels = char('estimate', 'st.err.', 't-stat.', 'p-value');
	if !silent
		printf('\n*********************************************************\n');
		printf('OLS estimation results\n');
		printf('Observations %d\n',n);
		printf('R-squared %f\n',rsq);
		printf('Sigma-squared %f\n',sigsq);
		p = 2 - 2*tcdf(abs(t), n - k);
		results = [b, seb, t, p];
		if regularvc
			printf('\nResults (Ordinary var-cov estimator)\n\n');
		else
			printf('\nResults (Het. consistent var-cov estimator)\n\n');
		end
		prettyprint(results, names, labels);
		printf('\n*********************************************************\n');
	end
end
