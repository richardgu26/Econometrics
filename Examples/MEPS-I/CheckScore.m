# This does simple Poisson estimation using the meps1996.data 

# The MEPS data

# Define dep and expl vbls
# 
# 	The dep. vbls, in corresponding column
# 	Office based doctor visits	    1
# 	Outpatient doctor visits		2
# 	Emergency room visits			3
# 	Inpatient visits				4
# 	Dental visits					5
# 	Prescriptions					6
# 	

1;


which_dep = 1;
if (which_dep == 1) printf("\nOBDV\n"); endif
if (which_dep == 2) printf("\nOPV");endif
if (which_dep == 3) printf("\nIPV");endif
if (which_dep == 4) printf("\nERV");endif
if (which_dep == 5) printf("\nDV");endif
if (which_dep == 6) printf("\nPRESCR");endif

load meps1996.data;
y = data(:,which_dep);
x = data(:,7:12);
n = rows(x);
x = [ones(n,1) x];

names = char("constant","pub. ins.","priv. ins.", "sex", "age","edu","inc");
model = "Poisson";
data = [y x];
theta = zeros(columns(x),1);
score = numgradient("Poisson", {theta, data}); 
printf("numeric gradient using unscaled data\n");
disp(mean(score));

x = scale_data(x);
data = [y x];
theta = zeros(columns(x),1);
score = numgradient("Poisson", {theta, data}); 
printf("numeric gradient using scaled data\n");
disp(mean(score));

