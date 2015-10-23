# Estimates the basic Nerlove Cobb-Douglas model
1;
function m = nm(theta, data)
	y = data(:,1);
	x = data(:,2:columns(data));
	e = y - x*theta;
	m = diag(e)*x;
endfunction

function m = nm2(theta, data)
	y = data(:,1);
	x = data(:,2:columns(data));
	e = y - x*theta;
	x = [x st_norm(x(:,2:5).^2)];
	m = diag(e)*x;
endfunction


load nerlove.data;

data = data(:,2:6);
data = log(data);
n = rows(data);
y = data(:,1);
x = data(:,2:5);
x = [ones(n,1), x];
%x = [ones(n,1), x, x(:,1).^2];

names = char("constant", "output","labor", "fuel", "capital","output2");
%names = char("constant", "output","labor", "fuel", "capital","output2");
data = [ y x ];

# estimate efficient weight matrix
b = mc_ols(y,x, "", 1);
e = y - x*b;
m = x .*e;
momentcov = m'*m/n;
weight = inv(momentcov);

# note that weights don't matter for results - exactly identified
gmm_results(b, data, 1, "nm", "", names, "Nerlove model estimated using GMM - Identity weight matrix", "", "", 0, momentcov);
gmm_results(b, data, weight, "nm", "", names, "Nerlove model estimated using GMM - Efficient weight matrix", "", "", 0, momentcov);

# now try an overidentified estimator
m = nm2(b, data);
momentcov = m'*m/n;
weight = inv(momentcov);
gmm_results(b, data, weight, "nm2", "", names, "Nerlove model estimated using GMM - Efficient weight matrix", "", "", 0, momentcov);


# OLS result for comparison
printf("OLS results for comparison\n");
mc_ols(y,x,names, 0, 1);
