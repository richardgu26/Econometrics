% creates moving average of lags of series. moving.m creates
% moving averages centered at period t. This does it centering
% so t-1 is the most recent observation used, going back l additionsl
% lags. For example, if l is 6, the return is a
% moving average of x(t-6) .. x(t-1)
function m = malags(x,l);
	m = lags(x,l);
	m = sum(m,2)/l;
endfunction	
