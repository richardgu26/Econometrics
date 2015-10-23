# this program illustrates sample selection bias associated with
# dropping observations for which the dep. vbl. is <0. The resulting figure is
# saved in sampsel.ps
clear all;
n = 100;
sig = 5;
x = 10*rand(n,1);
y = x+sig*randn(n,1);

z=y>0;
yin=y(z);  # this drops the rows for which y < 0
xin=x(z);

xx=[ones(size(xin)) xin];
b=inv(xx'*xx)*xx'*yin;


xx=(0:0.05:10);
yy=xx;
yhat=b(1)+b(2)*xx;

plot(x,y,'o;Data;',xx,yy,'--;True Line;',xx,yhat,'-;Fitted Line;');
print("sampsel.eps", "-depsc2");


