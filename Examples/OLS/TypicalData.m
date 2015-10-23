# this generates data according to classical model
# and shows the OLS fit

n = 20;
x = 1:n;
x = x';

x = [ones(n,1) x];
beta = [10; -1];
e = 6*randn(n,1);
true = x*beta;
y = true + e;

# the OLS coefficients
b = inverse(x'*x)*x'*y;

# Plot the fitted line
yhat = x*b;
x = x(:,2);
#title("Example OLS fit");
xlabel("X");
plot(x,y,"@;data;",x,true,"--;true regression line;");
pause(10);
%print("TypicalData.svg","-dsvg"); 
