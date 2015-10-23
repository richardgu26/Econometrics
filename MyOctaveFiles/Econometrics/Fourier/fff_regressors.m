# Sets up the regressor matrix for Gallant's Fourier form,
# omitting the quadratic terms.
#
# The J parameter has been incorporated into the kalpha
# argument. kalpha is a matrix of multi-indices, each
# multi-index is a row of the matrix. Since J directly
# multiplies the multi-indices, they are not required
# to be elementary: they may be positive integer multiples
# of elementary multi-indices.
#
# usage: z = fff_regressors(x, kalpha)

function z = fff_regressors(x, kalpha)

	minx = min(x);
	n = rows(x);
	x = x - repmat(minx, n, 1);
	maxx = max(x);
	x = x ./ repmat(maxx, n, 1);
	x = (2*pi-eps)*x;

	z = [ones(rows(x),1) x x.^2];

	# loop through multi-indices
	for i = 1:rows(kalpha)
		k = kalpha(i,:)';
		zz = x*k;
		zz = [2*cos(zz)  2*sin(zz)];
		z = [z zz];
	endfor

endfunction


