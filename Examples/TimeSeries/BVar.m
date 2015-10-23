% VAR with Minnesota priors
clc;

n = 30;
G = 3;
sig = 2*rand(3,3)-1;
sig = sig'*sig;
P = chol(sig);

% true 
A1 = [0.6 0.4 -0.3;
      -0.6 0.5 0.3;
      -0.5 0.3 0.7];
A2 = [0.1 0.1 -0.1;
      -0.1 -0.2 0.2;
      0.2 0 -0.3];
y1 = zeros(G,1);
y2 = zeros(G,1);
n = 100;
ys = zeros(n,G);
burnin = 100;
% generate the data
for t = 1:n+burnin
	y = A1*y1 + A2*y2 + P'*randn(G,1);
	y2 = y1;
	y1 = y;
	if (t> burnin) ys(t-burnin,:) = y'; endif
end

ylag = lags(ys,2);
ylag = ylag(3:end,:);
ys = ys(3:end,:);
b = ols(ys,ylag);
bols = b;
printf("OLS\n");
printf("matrix norm of A1-A1hat\n");
norm(A1-b(1:3,:)')
printf("matrix norm of A2-A2hat\n");
norm(A2-b(4:6,:)')

% VAR(2) with Minnesota priors
maxlag = 2;
g = columns(ys);
data = [ys lags(ys, maxlag)]; % add lags
data = data(maxlag+1:end,:); % drop rows with missing
y = data(:,1:g);
x = data(:,g+1:end);
n = rows(y);
% now Minnesota priors
PP = 1;   
y = [y; PP*eye(g)];     % first the restriction on first lag coefficients
if maxlag > 1           % now the higher order lags restricted to zero
        for j=2:maxlag
                y = [y; zeros(g,g)];
        endfor
endif
x = [x; PP*eye(g*maxlag)];
b = ols(y,x);
printf("BVAR\n");
printf("matrix norm of A1-A1hat\n");
norm(A1-b(1:3,:)')
printf("matrix norm A2-A2hat\n");
norm(A2-b(4:6,:)')


