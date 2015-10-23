# returns a stack of aux stats, each one in a row
# used for CUE-II (S cols. in randdraws)
function Z = aux_stat(theta, randdraws1, randdraws2)
	S = columns(randdraws1);
	Z = zeros(S, 6);
	for s = 1:S
		data = dgp(theta, randdraws1(:,s), randdraws2(:,s));
		b = data(:,1);
		% bound bid for numeric stability
		b = (b>0.01).*b + (b<0.01)*0.01;
		y = log(b); % use log of bid
		x = [ones(rows(y),1) data(:,2)];
		bhat = x\y;
		e = y - x*bhat;
		sig = log(e'*e/(rows(y)-2));
		m1 = mean(log(b));
		m2 = std(log(b));
		m3 = mean((log(b)-m1).^3);
		Z(s,:) = [bhat; sig; m1; m2; m3]';
	endfor
endfunction
