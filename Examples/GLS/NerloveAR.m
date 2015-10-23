# Estimates the basic Nerlove Cobb-Douglas model,
# and checks for autocorrelated errors by
# plotted residuals and using Breusch-Godfrey test

load nerlove.data;

data = data(:,2:6);
data = log(data);
n = rows(data);
y = data(:,1);
x = data(:,2:5);
x = [ones(n,1), x];

# get residuals
[junk, junk, e] = ols(y,x);

# quadratic fit to residuals
output = x(:,2);
xx = [ones(n,1) output output .^2];
[b, junk, junk] = ols(e, xx);
efit = xx*b;

plot(output, e, "o;Residuals;", output, efit, "-;Quadratic fit to Residuals;");
print("NerloveAR.eps", "-depsc2");

# do the test using 3 lags
BreuschGodfreyTest(e, 3, x);
