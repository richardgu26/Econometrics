# this generates data from the simple supply-demand model
# q = a1 + a2*p + a3*m + e1
# p = b1 + b2*p +        e2

# number of obsn
n = 100;
# model parameters
a1 = 100;
a2 = -1;
a3 = 1;
b1 = 20;
b2 = 1;
# exog var income
m = randn(n,5);
m = sum(m.*m,2); # chi square df=5
# structural shocks
e1 = 4*randn(n,1);
e2 = 4*randn(n,1);
# rf shocks
v1 = (b2*e1-a2*e2)/(b2-a2);
v2 = (e1-e2)/(b2-a2);
# rf coefs
pi11 = (b2*a1-a2*b1)/(b2-a2);
pi21 = b2*a3/(b2-a2);
pi12 = (a1-b1)/(b2-a2);
pi22 = a3/(b2-a2);
# sf variables
q = pi11 + pi21*m + v1;
p = pi12 + pi22*m + v2;
data = [q p m];
save -ascii data.txt data;


