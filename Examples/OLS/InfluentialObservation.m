# this illustrates effect and detection of influential observations

n = 20;
x = (1:n-1)/(n-1);
x = x';
x = [x; 3]; # the last observation is an outlying value of x

x = [ones(n,1) x];
P = x*inverse(x'*x)*x';

beta = [10; -1];
e = 2*randn(n,1);
y = x*beta + e;


# The fit
yhat = P*y;

# calculate leverage and influence
leverage = diag(P);
e = y - yhat;
influence = (leverage ./ (1-leverage)) .* e;

xlabel("X");
x = x(:,2);
plot(x, y, "o;Data points;", x, yhat, "-;fitted;", x, leverage, "--;Leverage;", x, influence, "--;Influence;");
axis([0,3.5]);
data = [y x];
save -ascii influencedata data;
print("InfluentialObservation.svg", "-dsvg");
