# computes autocorrelation matrices up to order p,
# and returns them vec'ed out from 1 to p
function rhos = autocors(x,p)
	n = rows(x);
	k = columns(x);
	for i = 1:p
		xx = x(1:n-p,:);
		xxx = x(p+1:n,:);
		rho = corr([xxx xx]);
		rho = rho(1:k,k+1:2*k);
		if i == 1; rhos = vec(rho);
		else rhos = [rhos; vec(rho)];
		endif
	endfor	
endfunction	
