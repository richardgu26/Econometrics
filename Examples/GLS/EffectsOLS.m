% this shows effects of het. on the OLS estimator
het = true;  % set this to true or false
n = 50; % sample size

% used by montecarlo.m to generate replications
function contrib = wrapper(args)
  het = args{1};
  n = args{2};
  x = randn(n,1);
  if het
    e = x .* randn(n,1); % note the heteroscedasticity here
  else
    e = randn(n,1); % note the heteroscedasticity here
  endif  
  y = x + e; % true coefs are zero for const, and 1 for slope
  x = [ones(n,1) x];
  [b, varb] = mc_ols(y,x, "", true, true);
  H0 = [0; 1]; % coefficient values under true null hypothesis
  t = (b - H0) ./ sqrt(diag(varb));
  contrib = t';
endfunction


args = {het, n};
reps = 1000;
n_pooled = 200;
system('rm data.effects'); % don't append, start from scratch
montecarlo('wrapper', args, reps, 'data.effects', n_pooled);

load data.effects;
t = abs(data(:,3:4));
crit_val = tinv(0.95, n-2); % 5% of prob. is to the R of this value
test = (t > crit_val) + (t < -crit_val); % now it's a 10% signif. level test
if het
  printf("data is heteroscedastic\n");
else
  printf("data is homoscedastic\n");
endif
printf('rejection frequency of nominal 10 percent test\n');
printf('intercept: %f\n', sum(test(:,1)/reps));
printf('slope: %f\n', sum(test(:,2)/reps));

bar([0 1], sum(test)/reps);
xlabel('0 is intercept, 1 is slope');
%print('EffectsOLS.png', '-dpng');

