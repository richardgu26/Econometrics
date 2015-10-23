function result = eemult_mv(m,v)

	if !ismatrix(m)
		error("eemult_mv: first arg must be a matrix");
	endif

	if !isvector(v)
		error("eemult_mv: second arg must be a vector");
	endif

	[rm, cm] = size(m);
	[rv, cv] = size(v);
	
	if (rm == rv)
		v = kron(v, ones(1,cm));
		result = m .* v;
	elseif (cm == cv)
		v = kron(v, ones(rm, 1));	
		result = m .* v;
	else
		error("eemult_mv: dimension of vector must match one of the dimensions of the matrix");	
	endif
endfunction




 
