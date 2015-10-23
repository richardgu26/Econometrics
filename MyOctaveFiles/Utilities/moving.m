%% Copyright (c) 2011, Aslak Grinsted
%% All rights reserved.
%%
%% Redistribution and use in source and binary forms, with or without
%% modification, are permitted provided that the following conditions are
%% met:
%%
%% * Redistributions of source code must retain the above copyright
%% notice, this list of conditions and the following disclaimer.
%% * Redistributions in binary form must reproduce the above copyright
%% notice, this list of conditions and the following disclaimer in
%% the documentation and/or other materials provided with the distribution
%%
%% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
%% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
%% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
%% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
%% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
%% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
%% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
%% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
%% POSSIBILITY OF SUCH DAMAGE.



function [y]=moving(x,m,fun)
%MOVING will compute moving averages of order n (best taken as odd)
%
%Usage: y=moving(x,n[,fun])
%where x 	is the input vector (or matrix) to be smoothed.
%      m 	is number of points to average over (best odd, but even works)
%      y 	is output vector of same length as x
%      fun  (optional) is a custom function rather than moving averages
%
% Note:if x is a matrix then the smoothing will be done 'vertically'.
%
%
% Example:
%
% x=randn(300,1);
% plot(x,'g.');
% hold on;
% plot(moving(x,7),'k');
% plot(moving(x,7,'median'),'r');
% plot(moving(x,7,@(x)max(x)),'b');
% legend('x','7pt moving mean','7pt moving median','7pt moving max','location','best')
%
% optimized Aslak Grinsted jan2004
% enhanced Aslak Grinsted Apr2007


if m==1
    y=x;
    return
end
if size(x,1)==1
    x=x';
end

if nargin<3
    fun=[];
elseif ischar(fun)
    fun=eval(['@(x)' fun '(x)']);
end

if isempty(fun)

    f=zeros(m,1)+1/m;
    n=size(x,1);
    isodd=bitand(m,1);
    m2=floor(m/2);


    if (size(x,2)==1)
        y=filter(f,1,x);
        y=y([zeros(1,m2-1+isodd)+m,m:n,zeros(1,m2)+n]);
    else
        y=filter2(f,x);
        y(1:(m2-~isodd),:)=y(m2+isodd+zeros(m2-~isodd,1),:);
        y((n-m2+1):end,:)=y(n-m2+zeros(m2,1),:);
    end

else
    y=zeros(size(x));
    sx=size(x,2);
    x=[nan(floor(m*.5),sx);x;nan(floor(m*.5),sx)];
    m1=m-1;
    for ii=1:size(y,1);
        y(ii,:)=fun(x(ii+(0:m1),:));
    end

end

return
