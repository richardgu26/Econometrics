%% Copyright (C) 2013 Michael Creel
%%
%% This file is free software; you can redistribute it and/or
%% modify it under the terms of the GNU General Public
%% License as published by the Free Software Foundation;
%% either version 3 of the License, or (at your option) any
%% later version.
%%
%% The code is distributed in the hope that it will be useful,
%% but WITHOUT ANY WARRANTY; without even the implied
%% warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
%% PURPOSE.  See the GNU General Public License for more
%% details.
%%
%% You should have received a copy of the GNU General Public
%% License along with Octave; see the file COPYING.  If not,
%% see <http://www.gnu.org/licenses/>.

% This example code shows how SBIL can be used to estimate the auction
% model in the paper "Indirect Likelihood Inference", by Michael Creel
% and Dennis Kristensen. This is a simple version that is for
% illustrative purposes. The actual code used to get the results for
% the paper implements parallelization, to allow work with large
% numbers of simulations. The intention of this code is to clearly
% show the ideas, without worrying about computational performance.
% To obtain results similar to those in the paper, increase S.
%
% usage: execute Auction_SBIL_Simple from the Matlab or Octave prompt.

function Auction_SBIL_Simple()
    S = 5e4;     	% number of draws from prior
    mc_reps = 1000;    	% number of Monte Carlo replication
    n = 80;		% sample size

    % initial parameter bounds: max mean value is 150, whatever the quality
    lb_param_ub = [
    -1     0.5  3   	% beta0
    0.0    0.5 	2  	% beta1
    ];

    theta0 = lb_param_ub(:,2); % true parameters
    lb = lb_param_ub(:,1);	% lower bound of prior
    ub = lb_param_ub(:,3);	% upper boudn of prior

    % generate the Monte Carlo replications at the design point
    Z_design = zeros(mc_reps, 6);
    for s = 1:mc_reps
        data = dgp(theta0, n);
    	Z_design(s,:) = auxstat(data);
    end

    % generate the S replications drawn from prior
    Z_paramspace = zeros(S, 6);
    theta_S = zeros(S, 2);
    for s = 1:S
	if(mod(s,1000)==0)
		fprintf('%d out of %d simulations done\n', s, S);
	end	
        theta_s = (ub-lb).*rand(2,1) + lb;
        theta_S(s,:) = theta_s';
        data = dgp(theta_s, n);
        Z_paramspace(s,:) = auxstat(data);
    end


    % number of neighbors to use
    k = floor(1.5*S^0.25);

    % find the indices of nearest neighbors
    idx = knnsearch(Z_design, Z_paramspace, k);

    % compute the posterior means
    for i = 1:size(Z_design,1)
        thetahat1(i,:) = mean(theta_S(idx(i,:),1));
        thetahat2(i,:) = mean(theta_S(idx(i,:),2));
    end

    % display basic results
    fprintf('SBIL estimation results\n');
    fprintf('Sample size: %d\n', n);
    fprintf('Number of simulated auxiliary statistics: %d\n', S);
    fprintf('Do not expect good results unless S > 10^5\n');

    fprintf('true parameter value theta1: %f\n', theta0(1,:));
    fprintf('true parameter value theta2: %f\n', theta0(2,:));
    fprintf('posterior mean theta1: %f\n', mean(thetahat1));
    fprintf('posterior mean theta2: %f\n', mean(thetahat2));
    fprintf('RMSE theta1: %f\n', sqrt(mean((thetahat1-theta0(1,:)).^2)));
    fprintf('RMSE theta2: %f\n', sqrt(mean((thetahat2-theta0(2,:)).^2)));
    hist(thetahat1, 50);
    title('thetahat1: true value is 0.5');
    figure;
    hist(thetahat2, 50);
    title('thetahat2: true value is 0.5');
end


function data = dgp(theta, n)
    % the model
    theta1 = theta(1,:);
    theta2 = theta(2,:);
    N = 6;
    % quality of good
    x = rand(n,1);
    % valuations drawn from exponetial mean phi
    phi = exp(theta1 + theta2*x);
    phi = repmat(phi,1,N); % there are N bidders with common valuation
    v = exprnd(phi);
    phi = phi(:,1);
    % highest valuation
    v = max(v')';
    % get winning bid
    z = v./phi;
    D = exp(-5*z).*(60*exp(5*z) + 300*phi .* exp(4*z) - 300*phi .* exp(3*z) ...
    + 200*phi .* exp(2*z) - 75*phi .* exp(z) + 12*phi)/60 - 137*phi/60;
    b = v - D ./ ((1 - exp(-v./phi)).^(N-1));
    % the data are the qualities and the winning bids
    data = [x b];
end


function Z = auxstat(data)
    x = data(:,1);
    b = data(:,2);
    % top bound bid to control large outliers
    b = b.*(b>0);
    % the aux stat
    b = 0.01.*(b<0.01) + b.*(b>0.01);
    y = log(b);
    z = [ones(size(x,1),1) x];
    bhat = z\y;
    e = y-z*bhat;
    s = log(e'*e/(size(y,1)-2));
    m1 = mean(log(b));
    m2 = std(log(b));
    m3 = mean((log(b)-m1).^3);
    Z = [bhat' s m1 m2 m3];
end



% the K nearest neighbors search used in this simple example uses
% the following code, which works fine on Matlab, which has JIT
% compiling. On Octave, it is slower. For large query and target
% sets, this code can become slow, and the ANN code is recommended.
% This is for illustrative purposes only.


% Copyright (c) 2009, Yi Cao
% All rights reserved.

% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:

%	* Redistributions of source code must retain the above copyright
%	notice, this list of conditions and the following disclaimer.
%	* Redistributions in binary form must reproduce the above copyright
%	notice, this list of conditions and the following disclaimer in
%	the documentation and/or other materials provided with the distribution
%
%	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%	AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%	IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
%	ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
%	LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
%	CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
%	SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
%	INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
%	CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
%	ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
%	POSSIBILITY OF SUCH DAMAGE.

function [idx,D]=knnsearch(varargin)
% KNNSEARCH   Linear k-nearest neighbor (KNN) search
% IDX = knnsearch(Q,R,K) searches the reference data set R (n x d array
% representing n points in a d-dimensional space) to find the k-nearest
% neighbors of each query point represented by eahc row of Q (m x d array).
% The results are stored in the (m x K) index array, IDX.
%
% IDX = knnsearch(Q,R) takes the default value K=1.
%
% IDX = knnsearch(Q) or IDX = knnsearch(Q,[],K) does the search for R = Q.
%
% Rationality
% Linear KNN search is the simplest appraoch of KNN. The search is based on
% calculation of all distances. Therefore, it is normally believed only
% suitable for small data sets. However, other advanced approaches, such as
% kd-tree and delaunary become inefficient when d is large comparing to the
% number of data points. On the other hand, the linear search in MATLAB is
% relatively insensitive to d due to the vectorization. In  this code, the
% efficiency of linear search is further improved by using the JIT
% aceeleration of MATLAB. Numerical example shows that its performance is
% comparable with kd-tree algorithm in mex.
%
% See also, kdtree, nnsearch, delaunary, dsearch

% By Yi Cao at Cranfield University on 25 March 2008

% Example 1: small data sets
%{
R=randn(100,2);
Q=randn(3,2);
idx=knnsearch(Q,R);
plot(R(:,1),R(:,2),'b.',Q(:,1),Q(:,2),'ro',R(idx,1),R(idx,2),'gx');
%}

% Example 2: ten nearest points to [0 0]
%{
R=rand(100,2);
Q=[0 0];
K=10;
idx=knnsearch(Q,R,10);
r=max(sqrt(sum(R(idx,:).^2,2)));
theta=0:0.01:pi/2;
x=r*cos(theta);
y=r*sin(theta);
plot(R(:,1),R(:,2),'b.',Q(:,1),Q(:,2),'co',R(idx,1),R(idx,2),'gx',x,y,'r-','linewidth',2);
%}

% Example 3: cputime comparion with delaunay+dsearch I, a few to look up
%{
R=randn(10000,4);
Q=randn(500,4);
t0=cputime;
idx=knnsearch(Q,R);
t1=cputime;
T=delaunayn(R);
idx1=dsearchn(R,T,Q);
t2=cputime;
fprintf('Are both indices the same? %d\n',isequal(idx,idx1));
fprintf('CPU time for knnsearch = %g\n',t1-t0);
fprintf('CPU time for delaunay  = %g\n',t2-t1);
%}
% Example 4: cputime comparion with delaunay+dsearch II, lots to look up
%{
Q=randn(10000,4);
R=randn(500,4);
t0=cputime;
idx=knnsearch(Q,R);
t1=cputime;
T=delaunayn(R);
idx1=dsearchn(R,T,Q);
t2=cputime;
fprintf('Are both indices the same? %d\n',isequal(idx,idx1));
fprintf('CPU time for knnsearch = %g\n',t1-t0);
fprintf('CPU time for delaunay  = %g\n',t2-t1);
%}
% Example 5: cputime comparion with kd-tree by Steven Michael (mex file)
% <a href="http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=7030&objectType=file">kd-tree by Steven Michael</a>
%{
Q=randn(10000,10);
R=randn(500,10);
t0=cputime;
idx=knnsearch(Q,R);
t1=cputime;
tree=kdtree(R);
idx1=kdtree_closestpoint(tree,Q);
t2=cputime;
fprintf('Are both indices the same? %d\n',isequal(idx,idx1));
fprintf('CPU time for knnsearch = %g\n',t1-t0);
fprintf('CPU time for delaunay  = %g\n',t2-t1);
%}

% Check inputs
[Q,R,K,fident] = parseinputs(varargin{:});

% Check outputs
error(nargoutchk(0,2,nargout));

% C2 = sum(C.*C,2)';
[N,M] = size(Q);
L=size(R,1);
idx = zeros(N,K);
D = idx;

if K==1
% Loop for each query point
for k=1:N
d=zeros(L,1);
for t=1:M
d=d+(R(:,t)-Q(k,t)).^2;
end
if fident
d(k)=inf;
end
[D(k),idx(k)]=min(d);
end
else
for k=1:N
d=zeros(L,1);
for t=1:M
d=d+(R(:,t)-Q(k,t)).^2;
end
if fident
d(k)=inf;
end
[s,t]=sort(d);
idx(k,:)=t(1:K);
D(k,:)=s(1:K);
end
end
if nargout>1
D=sqrt(D);
end

end

function [Q,R,K,fident] = parseinputs(varargin)
% Check input and output
error(nargchk(1,3,nargin));

Q=varargin{1};

if nargin<2
R=Q;
fident = true;
else
fident = false;
R=varargin{2};
end

if isempty(R)
fident = true;
R=Q;
end

if ~fident
fident = isequal(Q,R);
end

if nargin<3
K=1;
else
K=varargin{3};
end

end
