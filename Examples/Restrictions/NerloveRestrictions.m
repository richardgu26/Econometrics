# Tests CRTS and HOD1 for the Nerlove Cobb-Douglas model

load nerlove.data;

data = data(:,2:6);
data = log(data);
n = rows(data);
y = data(:,1);
x = data(:,2:5);
x = [ones(n,1), x];

names = char("constant", "output","labor", "fuel", "capital");

# First HOD1
R = [0, 0, 1, 1, 1];
r = 1;
printf("\nImposing and testing HOD1\n");
mc_olsr(y, x, R, r, names);
TestStatistics(y, x, R, r);

# Now CRTS
R = [0, 1, 0, 0, 0];
r = 1;
printf("\nImposing and testing CRTS\n");
mc_olsr(y, x, R, r, names);
TestStatistics(y, x, R, r);
