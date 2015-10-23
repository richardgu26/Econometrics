# this generates data according to classical model
# and shows the OLS fit

n = 20;
x = 1:n;
x = x';

x = [ones(n,1) x];
beta = [12; -1.5];
true = x*beta;
e = 4*randn(n,1);
y = true + e;

# the OLS coefficients
b = inverse(x'*x)*x'*y;

# Plot the fitted line
yhat = x*b;
x = x(:,2);
#title("Example OLS fit");
xlabel("X");
plot(x, y, "*;data points;", x, yhat, "-;fitted line;", x, true, "-;true line;")
%print("OlsFit.svg", "-dsvg");


