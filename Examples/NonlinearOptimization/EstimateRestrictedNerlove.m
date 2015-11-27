% Nerlove model by OLS, but imposing that factor shares are in (0,1) and sum to 1
load nerlove.data
n = size(data,1);
y = log(data(:,2));
x = [ones(n,1) log(data(:,3:6))];

thetastart = [zeros(5,1)];

% unrestricted OLS using fminunc
[thetahat, ssr] = fminunc(@(theta) (y-x*theta)'*(y-x*theta), thetastart);
fprintf('the OLS estimates\n');
thetahat
ssr

% restricted OLS using sqp
% the box bounds
lb = [-Inf; 0; 0; 0; 0];
ub = [Inf; Inf;1;1;1];
% the equality restriction: bL + b_f + bK = 1
R = [0 0 1 1 1];
r = 1;       
[thetahat, ssr] = sqp(thetastart, @(theta) (y-x*theta)'*(y-x*theta), @(theta) R*theta-r, [], lb, ub); %Octave replacement for fmincon
fprintf('\nrestricted estimates\n');
thetahat
ssr
