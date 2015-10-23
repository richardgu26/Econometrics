// Gurmu MGF (mean restricted version only)
//  This MGF is used to define the Gurmu semiparametric count density.
//  We'll use psi internally as the entire parameter to follow Gurmu's notation
//  set order = zeros(n,1) for binary part of hurdle
//  set order = y for original SP or count part of hurdle SP

#include <octave/oct.h>


DEFUN_DLD(GurmuMGF, args, ,"Gurmu MGF, C++ version\nmgf = GurmuMGF(theta, x, order)\n\
where theta holds the parameters in the order\n\
theta = [beta;alpha;c]\n\
It is a assumed that alpha>0 has been imposed in calling program\n\
Also, this imposes the restriction that the mean of y is exp(x*beta)")

{
	Matrix psi(args(0).matrix_value());
	Matrix x(args(1).matrix_value());
	Matrix order(args(2).matrix_value()); // the vector y goes here, or a vector of zeros for GurmuBinary

	// parameters are beta, alpha, c.
	const int n = args(1).rows();
	const int p = args(0).rows();
	int i, j, k, l, m;
	double temp, alpha,term1, lambda, combo1, combo2, hj, hk;

 	k = args(1).columns();
	alpha = psi(k,0); // 0 < alpha already imposed in calling functions

	Matrix beta(k,1);
	Matrix c(p-k,1);
	Matrix theta;
	Matrix term2(n,1);
	Matrix term3(n,1);
	Matrix gurmu_moment(n,1);


	for (i=0;i<k;i++)
	{
		beta(i,0) = psi(i,0);
	}	
	
	c(0,0) = 1;
	for (i = 1;i < c.rows();i++)
	{
		c(i,0) = psi(k+i,0);
	}	

	theta = x*beta;
	for (i=0; i < n; i++)
	{
		temp = theta(i,0);
		theta(i) = exp(temp);
	}
	
	temp = 0; // calculate c'c
	for (i=0; i < c.rows(); i++)
	{
		temp = temp + c(i,0)*c(i,0);
	}
	
	term1 = exp(lgamma(alpha)) / temp;

	lambda = 0; 

	for (j=0; j < c.rows(); j++)
	{
		for (k=0; k < c.rows(); k++)
		{
			for (l=0; l <= j; l++)
			{
				for (m=0; m <= k; m++)
				{
					hj = lgamma(j + alpha) - lgamma(alpha) - lgamma(j + 1);
					hk = lgamma(k + alpha) - lgamma(alpha) - lgamma(k + 1);
					hj = exp(hj);
					hk = exp(hk);
					combo1 = lgamma(j + 1) - lgamma(l + 1) - lgamma(j - l + 1);
					combo1 = exp(combo1);
					combo2 = lgamma(k + 1) - lgamma(m + 1) - lgamma(k - m + 1);
					combo2 = exp(combo2);
					lambda = lambda + c(j,0) * c(k,0) * sqrt(hj*hk) \
						* combo1 * combo2 * exp(lgamma(alpha + l + m + 1) - lgamma(alpha + l) - lgamma(alpha + m)) \
					 	* pow(-1.0,-((double) l + (double) m));
				}
			}
		}
	}


	lambda = lambda * term1;  // last part of Gurmu eqn 17

	for (i=0; i < n; i++)
	{
		term2(i,0) = pow(1.0 + theta(i,0)/lambda, -alpha) * pow(lambda + theta(i,0), -order(i,0));
	}
	
	term3.fill(0.0);
	for (j=0; j < c.rows(); j++)
	{
		for (k=0; k < c.rows(); k++)
		{
			for (l=0; l <= j; l++)
			{
				for (m=0; m <= k; m++)
				{
					hj = lgamma(j + alpha) - lgamma(alpha) - lgamma(j + 1);
					hk = lgamma(k + alpha) - lgamma(alpha) - lgamma(k + 1);
					hj = exp(hj);
					hk = exp(hk);
					combo1 = lgamma(j+1) - lgamma(l+1) - lgamma(j - l + 1);
					combo1 = exp(combo1);
					combo2 = lgamma(k+1) - lgamma(m+1) - lgamma(k - m + 1);
					combo2 = exp(combo2);
					for (i=0; i<n; i++)
					{
						term3(i) = term3(i) + c(j) * c(k) * sqrt(hj*hk) \
							* combo1 * combo2 \
							* exp(lgamma(alpha + l + m + order(i)) - lgamma(alpha + l) - lgamma(alpha + m)) \
							* pow((-1.0 - theta(i)/lambda), -(double) l - (double) m);
					}	
				}
			}
		}
	}


	for (i=0; i<n; i++)
	{	
		gurmu_moment(i,0) = term1 * term2(i,0) * term3(i,0);
  	}
	

	return octave_value(gurmu_moment);
}
