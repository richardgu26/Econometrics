usematlab = false;

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
thetastart = [mean(y); var(y); 0.1];
[thetahat, logL] = fminunc(@(theta) -arch1(theta, data), thetastart);

logL = -logL; % re-convert

%%%%%%%%%%%%%%%%   Results   %%%%%%%%%%%%%%%%
fprintf('ARCH(1) results for NYSEWK data set (fminunc)\n');

fprintf('the estimate\n');
disp(thetahat);

fprintf('the logL value\n');
disp(logL);

BIC = -2*logL + size(thetahat,1)*log(size(y,1));
fprintf('BIC value\n');
disp(BIC);

%%
%%%%%%%%%%%%%%% Constrained maximization of logL  %%%%%%%%%%%
%  ARCH model needs parameter restrictions for stationariry and positive variance
%  fmincon can be used to impose them
lb = [-Inf; 0; 0];
ub = [Inf; Inf; 1];
if usematlab
  [thetahat, logL] = fmincon(@(theta) -arch1(theta, data), thetastart, [],[], [],[], lb, ub);
else
  [thetahat, logL] = sqp(thetahat, @(theta) -arch1(theta,data),[], [], lb, ub); %Octave replacement for fmincon
end

logL_UR = -logL; % re-convert

%%%%%%%%%%%%%%%%   Results   %%%%%%%%%%%%%%%%
fprintf('ARCH(1) results for NYSEWK data set (fmincon - impose stationarity) \n');

fprintf('the estimate\n');
disp(thetahat);

fprintf('the logL value\n');
disp(logL_UR);

%%
%%%%%%%% here's an example of a binding constraint: use the results to do LR test   
R = [0 1 1]; % a silly constraint: make last 2 coefficients add to 2.5
r = 2.5;
lb = [-Inf; 0; 0];
ub = [Inf; Inf; 1];
if usematlab
  [thetahat_r, logL] = fmincon(@(theta) -arch1(theta, data), thetastart, [],[],R,r, lb, ub);
else
  [thetahat_r, logL] = sqp(thetahat, @(theta) -arch1(theta,data),@(theta) R*theta-r, [], lb, ub); %Octave replacement for fmincon
end
logL_R = -logL; % re-convert

%%%%%%%%%%%%%%%%   Results   %%%%%%%%%%%%%%%%
fprintf('ARCH(1) results for NYSEWK data set (fmincon - impose silly restriction) \n');

fprintf('the estimate\n');
disp(thetahat_r);

fprintf('the logL value\n');
disp(logL_R);

fprintf('the LR test statistic\n');
disp(2*(logL_UR-logL_R));

