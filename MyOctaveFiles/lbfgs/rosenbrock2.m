function [f, df] = rosenbrock2(x);

% rosenbrock.m This function returns the function value, partial derivatives
% and Hessian of the (general dimension) rosenbrock function, given by:
%
%       f(x) = sum_{i=1:D-1} 100*(x(i+1) - x(i)^2)^2 + (1-x(i))^2 
%
% where D is the dimension of x. The true minimum is 0 at x = (1 1 ... 1).
%
% Carl Edward Rasmussen, 2001-07-21.

D = length(x);
f = sum(100*(x(2:D)-x(1:D-1).^2).^2 + (1-x(1:D-1)).^2);

if nargout > 1
 df = numgradient("rosenbrock2",{x});
end

