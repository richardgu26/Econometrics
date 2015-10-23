# decriptive statistics useful for count data
# call with 2 args to see screen output, 
# or with only the data to get output returned to you

function [info, rlabels] = CountDstats(x, screenout);

test = x == 0;
perc_zero = sum(test)/rows(x);
info = [mean(x); std(x); (mean(x) ./ std(x) .^ 2); min(x); max(x); perc_zero]; 

rlabels = char("mean","st. dev.","mean/var", "min","max","\\%% zero");

if nargin == 2 prettyprint_r(info, rlabels); endif

endfunction
