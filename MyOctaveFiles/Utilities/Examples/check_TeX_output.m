Utilities;
 
a = rand(3,4);
b = rand(2,1);
b = postpad(b,rows(a),nan);
data = [a, b];
sterr = rand(size(data));
rname = char("r1","r2","r3");
cname = char("c1","c2","c3","c4","c5");

results2tex(data, sterr, "test.tex", cname, rname, 5, "r", 0);
