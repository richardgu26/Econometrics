# Estimates 5 separate models for the Nerlove Cobb-Douglas model,
# and does a Chow test for pooled coefficients

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
	big_x(startrow:endrow,startcol:endcol) ...
		= big_x(startrow:endrow,startcol:endcol) ...
		+ x(startrow:endrow,:);
endfor
x = big_x;

names = char("constant", "output","labor", "fuel", "capital");
names = [names; names; names; names; names]; # copy 5 times

printf("\nNerlove model: 5 separate regressions for different output levels\n");
b = mc_ols(y, x, names);
output = 1:5;
output = output*5+2-5;
output = b(output,:);
rts = 1 ./ output;

# gset term X11
xlabel("Output group");
group = 1:5;
group = group';
plot(group, rts,"-;RTS;")
print("rts.svg","-dsvg");


# Chow test
R = eye(5);
Z = zeros(5,5);
R = [
	R -R Z Z Z;
	R Z -R Z Z;
	R Z Z -R Z;
	R Z Z Z -R
	];
r = zeros(20,1);

printf("\nChow test: note that the restricted model\n");
printf("gives the same results as the original model\n");
mc_olsr(y, x, R, r, names);
TestStatistics(y, x, R, r);



