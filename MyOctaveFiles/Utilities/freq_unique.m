## freq_unique: counts frequencies of unique elements of a vector


function a = freq_unique(x)
	b = unique(x);
	n = rows(b);
	a = zeros(n,1);
	for i=1:n
		count = b(i) == x;
		a(i)=sum(count);
	endfor
	a = sort(a,1,'descend');
endfunction	
