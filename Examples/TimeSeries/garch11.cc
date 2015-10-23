// garch11 model, to illustrate using C++ for Octave bottlenecks
#include <oct.h>
#include <octave/Cell.h>

DEFUN_DLD(garch11, args, ,"garch11")
{
	// parameter of model
	ColumnVector model_params (args(0).column_vector_value());
	const double mu = model_params(0);
	const double omega = model_params(1);
	const double alpha = model_params(2);
	const double beta = model_params(3);

	// data
	const ColumnVector y (args(1).column_vector_value());

	const int n = y.rows();
	double h, hlag, e, elag, c, logdensity;
	const double pi = 3.1415926;
	int t;
	octave_value_list f_return;

	c = -log(sqrt(2.0*pi)); // constant part of log likelihood

	elag = 0.0;

	// initialize variance as variance of first 10 obs.
        hlag = 0.0;
        for (t=0; t<10;t++) hlag += (y(t)-mu)*(y(t)-mu);
        hlag = hlag/9.0;

	for (t=0; t < n; t++) { // loop over time
                e = y(t) - mu;
		h = omega + alpha*elag*elag + beta*hlag;
		logdensity += c -0.5*log(h) - 0.5*e*e/h;
		hlag = h;
		elag = e;
	}
	f_return(0) = logdensity;
	return f_return;
}
