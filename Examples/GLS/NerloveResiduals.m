# Nerlove Cobb-Douglas model with constant and output varying, input price coefs fixed

load nerlove.data;

data = data(:,2:6);
data = log(data);
n = rows(data);
y = data(:,1);
x = data(:,2:5);
x = [ones(n,1), x];
k = columns(x);

# create the block diagonal X matrix corresponding to separate coefficients
big_x = zeros(n,5*k);
for i=1:k
	startrow = (i-1)*29+1;
	endrow = i*29;
	startcol =(i-1)*k + 1;
	endcol = i*k;
	big_x(startrow:endrow,startcol:endcol) \
		= big_x(startrow:endrow,startcol:endcol) \
		+ x(startrow:endrow,:);
endfor
x = big_x;

names = char("constant", "output","labor", "fuel", "capital");
names = [names; names; names; names; names]; # copy 5 times

# Now let's try out the model with input prices restricted but constant and output varying
R = eye(5);
R = R(3:5,:);
Z = zeros(3,5);
R = [
	R -R Z Z Z;
	R Z -R Z Z;
	R Z Z -R Z;
	R Z Z Z -R
	];

r = zeros(12,1);

printf("Constant and output vary across the groups, input prices fixed\n");
mc_olsr(y, x, R, r, names);
TestStatistics(y, x, R, r);

PlotResiduals(y,x,1);
pause(10);

