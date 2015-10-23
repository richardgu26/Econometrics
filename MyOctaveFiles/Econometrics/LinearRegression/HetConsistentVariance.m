function varb = HetConsistentVariance(x, e)
	xx_inv = inv(x'*x);
	E = e .^2;
	varb = xx_inv * x'*eemult_mv(x, E) * xx_inv;
end
