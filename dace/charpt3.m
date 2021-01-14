function charpt3 
x0=-1:0.5:2;
y0=[-4.447 -0.452 0.551 0.048 -0.447 0.549 4.552];
n=3;%n is the order of regression
alph=polyfit(x0,y0,n);% regression funtion belongs to the software
y=polyval(alph,x0);% to solve the function value of x0
r=(y0-y)*(y0-y)';%mean square
x=-1:0.01:2;
y=polyval(alph,x);
plot(x,y,'k--');
xlabel('x');ylabel('y0*and polyfit.y--');
hold on;
plot(x0,y0,'*');
title('polynomial regression for distributed data');
grid on;
disp(['mean square£º',sprintf('%g\n',r)]);
disp(['parameter alph£º',sprintf('%g\t',alph)]);
c=3;
d=f(c);
d

end

