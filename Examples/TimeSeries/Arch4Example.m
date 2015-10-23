usematlab = false;
% Arch1Example results for start values
thetastart = [0.16; 3.17; 0.25; 0;0;0];

%%
% weekly close price of NSYE, data provided with GRETL
load nysewk.mat;
n = size(nysewk);
y = 100 * log( nysewk(2:n) ./ nysewk(1:n-1) );
data = y;
plot(y);


%%
%%%%%%%%%%%%%%% Constrained maximization of logL  %%%%%%%%%%%
%  ARCH model needs parameter restrictions for stationariry and positive variance
%  fmincon can be used to impose them
lb = [-Inf; 0; 0; 0; 0; 0];
ub = [Inf; Inf; 1;1;1;1];
R = [0 0 1 1 1 1]; % 
r = 1;
if usematlab
  [thetahat, logL] = fmincon(@(theta) -arch4(theta, data), thetastart, R,r, [],[], lb, ub);
else
  [thetahat, logL] = sqp(thetastart, @(theta) -arch4(theta,data),[], @(theta) r-R*theta, lb, ub); %Octave replacement for fmincon
end
logL_UR = -logL; % re-convert
i
%%%%%%%%%%%%%%%%   Results   %%%%%%%%%%%%%%%%
fprintf('ARCH(4) results for NYSEWK data set (fmincon - impose stationarity) \n');

fprintf('the estimate\n');
disp(thetahat);

fprintf('the logL value');
disp(logL_UR);


BIC = -2*logL + size(thetahat,1)*log(size(y,1));
fprintf('BIC value');
disp(BIC);
