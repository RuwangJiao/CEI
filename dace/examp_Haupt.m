a=[3 3];
b=[0 0;3.7 4];
XI=gridsamp(b,a);
Y1=cfun1(XI);
Y2=cfun2(XI);
global dmodel1;
global dmodel2;
theta = [10 10]; lob = [1e-1 1e-1]; upb = [20 20];
[dmodel1, perf] =dacefit(XI, Y1, @regpoly1, @corrgauss, theta, lob, upb);
[dmodel2, perf] =dacefit(XI, Y2, @regpoly1, @corrgauss, theta, lob, upb);
x=[3.0,2.6];
[x fval] = ga(@cbs,2);


% X=gridsamp([0 0;3.7 4],20);
% [YX,MES]=predictor(X,dmodel);
% 
% X1 = reshape(X(:,1),20,20); 
% X2 = reshape(X(:,2),20,20);
% YX = reshape(YX, size(X1));
% mesh(X1,X2,YX);figure(2),
% mesh(X1,X2,reshape(MES,size(X1)));

function cfun1=cfun1(x)
cfun1=-x(1).*sin(4*x(1))-1.1*x(2).*sin(2*x(2));
end
function cfun2=cfun2(x)
cfun2=x(1)+x(2)-3;
end
function cbs=cbs(x)
global dmodel1;
global dmodel2;
[cf1,cf1mse]=predictor(x,dmodel1);
[cf2,cf2mse]=predictor(x,dmodel2);
if cf1>=0 & cf2>=0
    cbs=(normpdf(cf1/cf1mse,0,1)+normpdf(cf2/cf2mse,0,1));
else
    cbs=0;
end
end