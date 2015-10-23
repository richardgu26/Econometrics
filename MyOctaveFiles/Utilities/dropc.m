## usage: 	b = dropc(m,c)
## inputs: 
##	m: a matrix
## 	c: a scalar
## outputs:
##	b: the matrix m, having dropped column c

function m = dropc(m,c) 
	k = columns(m);
	if c==1
		m = m(:,2:k);
	elseif c==k
		m = m(:,1:k-1);
	else
		m = [m(:,1:c-1) m(:,c+1:k)];
	endif
endfunction		 
	
