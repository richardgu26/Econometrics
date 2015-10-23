a = rand(2,1);
b = rand(3,3);
a = postpad(a,rows(b),nan);
data = [a, b];
sterr = rand(size(data));
rname = char("r1","r2","r3");
cname = char("c1","c2","c3","c4");

results2tex(data, sterr, "test.tex", cname, rname, 5, "r", 1);

