%%
% weekly close price of NSYE, data provided with GRETL
load nysewk.mat;
n = size(nysewk);
y = 100 * log( nysewk(2:n) ./ nysewk(1:n-1) );
data = y;
plot(y);

%%
%%%%%%%%%%%%%%% Unconstrained maximization of logL  %%%%%%%%%%%
% note that the objective has a minus sign in front, as fminunc
% minimizes, but we want to maximize the logL
thetastart = [mean(y); var(y); 0.1; 0.1];
[thetahat, logL] = fminunc(@(theta) -garch11(theta, data), thetastart);

logL = -logL; % re-convert

%%%%%%%%%%%%%%%%   Results   %%%%%%%%%%%%%%%%
fprintf('GARCH(1,1) results for NYSEWK data set (fminunc)\n');

fprintf('the estimate\n');
disp(thetahat);

fprintf('the logL value');
disp(logL);

BIC = -2*logL + size(thetahat,1)*log(size(y,1));
fprintf('BIC value');
disp(BIC);
