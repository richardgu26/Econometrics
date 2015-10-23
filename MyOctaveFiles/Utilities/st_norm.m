	# standardize and normalize data
	function [zz, m, s] = st_norm(z);
		n = rows(z);
		m = mean(z);
		s = std(z);
		zz = (z - repmat(m, n, 1)) * diag(1 ./ s);
	endfunction
