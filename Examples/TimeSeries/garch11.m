% computes GARCH(1,1) log likelihood
function logL = garch11(theta, data)
    mu = theta(1,:);
    omega = theta(2,:);
    alpha = theta(3,:);
    beta = theta(4,:);
    y = data(:,1);
    e = y - mu;
    n = size(e,1);
    h = zeros(n,1);
    
    % either of these next two are reasonable choices
    h(1,:) = var(y(1:10,:));
    %h(1,:) = var(y);
    for t = 2:n
        h(t,:) = omega + alpha*e(t-1,:)^2 + beta*h(t-1,:);
    end
    logL = -log(sqrt(2*pi)) -0.5*log(h) - 0.5*(e.^2)./h;
    logL = sum(logL);
