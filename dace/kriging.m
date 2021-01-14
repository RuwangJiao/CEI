function kriginghdmr
load data1
theta=[10 10];lob=[1e-1 1e-1];upb=[20 20];
[dmodel,perf]=...
dacefit(S,Y,@regpoly0,@corrgauss,theta,lob,upb)
X=gridsamp([0 0;100 100],40);
[YX MSE]=predictor(X,dmodel);
X1=reshape(X(:,1),40,40);X2=reshape(X(:,2),40,40);
YX=reshape(YX,size(X1));
figure(1),mesh(X1,X2,YX)
hold on,
plot3(S(:,1),S(:,2),Y,'.k','MarkerSize',10)
hold off


end

