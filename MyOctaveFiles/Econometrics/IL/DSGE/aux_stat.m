function Z = aux_stat(data)
	# variables coming in are y c n r w (output cons hours intrate wages)
	if (all(std(data) > 0))

		% bound data
		test = data > 1000;
		data = test*1000+ (1-test).*data;

		logdata = log(data);
		lagdata = lags(data,1);
		laglogdata = lags(logdata,1);
		n = rows(data);
		data = data(2:n,:);
		logdata = logdata(2:n,:);
		lagdata = lagdata(2:n,:);
		laglogdata = laglogdata(2:n,:);

		checkme = [data lagdata logdata laglogdata];

		% check for bad inputs
		if (sum(any(isnan(checkme)))|| sum(any(isinf(checkme))))
			Z = -1000*ones(7 + 5 + 15 + 25,1);
		else

			output = data(:,1);	% output
			cons = data(:,2);	% consumption
			hours = data(:,3);	% hours
			intrate = data(:,4);	% intrate
			wages = data(:,5);	% wages

			logoutput = logdata(:,1);	% output
			logcons = logdata(:,2);		% consumption
			loghours = logdata(:,3);	% hours
			logintrate = logdata(:,4);	% intrate
			logwages = logdata(:,5);	% wages

			lagoutput = lagdata(:,1);	% output
			lagcons = lagdata(:,2);		% consumption
			laghours = lagdata(:,3);	% hours
			lagintrate = lagdata(:,4);	% intrate
			lagwages = lagdata(:,5);	% wages

			nobs = rows(data);

			% alpha, rho1, sig1 (production function)
			% use a MA in investment as proxy for capital
			investment = output - cons;
			capitalproxy = investment(1:nobs-11,:) +...
				investment(2:nobs-10,:)+...
 				investment(3:nobs-9,:)+...
				investment(4:nobs-8,:)+...
				investment(5:nobs-7,:)+...
				investment(6:nobs-6,:)+...
				investment(7:nobs-5,:)+...
				investment(8:nobs-4,:)+...
				investment(9:nobs-3,:)+...
				investment(10:nobs-2,:);
 			logcapitalproxy = log(capitalproxy);
			x = [logcapitalproxy loghours(12:nobs,:)];
			x = [ones(rows(x),1) x];
			y = logoutput(12:nobs,:);
			w = laglogdata(12:nobs,:);
			w = [ones(rows(w),1) w];
			W = w*inverse(w'*w)*w';
			xhat = W*x;
			b = xhat\y;
			%b = x\y;
			e = y-x*b;
			junk = [e lag(e,1)];
			junk = junk(2:rows(junk),:);
			y = junk(:,1);
			x = junk(:,2);
			rho1 = corr(x,y);
			e = y-x*rho1;
			sig1 = e'*e/nobs;
			Z = [1-b(3,:); rho1; sig1];


			% gam, rho2, sig2 (MRS=wage)
			x = [ones(nobs,1) logcons];
			y = logwages;
			w = laglogdata;
			w = [ones(rows(w),1) w];
			W = w*inverse(w'*w)*w';
			xhat = W*x;
			b = xhat\y;
			%b = x\y;
			e = y-x*b;
			junk = [e lag(e,1)];
			junk = junk(2:rows(junk),:);
			y = junk(:,1);
			x = junk(:,2);
			rho2 = corr(y,x);
			e = y-x*rho2;
			sig2 = e'*e/nobs;
			Z = [Z; b; rho2; sig2];

			# standard devs. and correlations
			data = data(:,1:5);
			s = std(data);
			s2  = corr(data);
			s2 = s2 - diag(diag(s2)); % remove 1s on main diag
			s2 = diag(s) + s2;
			s2 = vech(s2);

			% autocorrelations
			m = mean(data);
			temp = (data - m)./s;
			temp2 = (lagdata - m)./s;
			ac = mean((temp).*temp2);
			for i = 1:4
				temp = shift(temp,i,2);
				ac2 = mean((temp).*temp2);
				ac = [ac ac2];
			endfor
			ac = ac';
			Z = [Z; m'; s2; ac];
			if sum(any(isnan(Z)) || any(isinf(Z)))
				Z = -1000*ones(7 + 5 + 15 + 25,1);
			endif
		endif
	else
		Z = -1000*ones(7 + 5 + 15 + 25,1);
	endif
endfunction


