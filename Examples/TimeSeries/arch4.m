% computes ARCH(4) log likelihood
function logL = arch4(theta, data)
    mu = theta(1,:);
    omega = theta(2,:);
    alpha1 = theta(3,:);
    alpha2 = theta(4,:);
    alpha3 = theta(5,:);
    alpha4 = theta(6,:);
    y = data(:,1);
    e = y - mu;
    n = size(e,1);
    h = zeros(n,1);
    
    % either of these next two are reasonable choices
    % h(1:4,:) = var(y(1:10,:));
    h(1:4,:) = var(y);
    for t = 5:n
        h(t,:) = omega + alpha1*e(t-1,:)^2 + alpha2*e(t-2,:)^2+ alpha3*e(t-3,:)^2+ alpha4*e(t-4,:)^2;
    end
    logL = -log(sqrt(2*pi)) -0.5*log(h) - 0.5*(e.^2)./h;
    logL = sum(logL);