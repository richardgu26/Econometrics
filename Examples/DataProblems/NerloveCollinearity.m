# Check the Nerlove data for collinearity
load nerlove.data;

data = data(:,3:6); # these are the regressors
data = log(data);
n = rows(data);

# Indices of remaining regressors, after selecting each in turn
j = [2 3 4; 1 3 4; 1 2 4; 1 2 3];
title = char("output", "labor","capital","fuel");
for i = 1:4
	printf("\n");
	disp(title(i,:));
	y = data(:,i); # select the artificial dep. var.
	k = j(i,:); # select the other regressors
	names = ["constant"; title(k,:)];
	x = [ones(n,1) data(:,k)];
	mc_ols(y, x, names);
endfor	
