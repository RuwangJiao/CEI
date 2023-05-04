function [ 	N_1 ] = newsecond4(m, N_1 )
%NEWSECOND2 Summary of this function goes here
%   Detailed explanation goes here
C=lhsamp(m,2)
C1=C(:,1);
C2=C(:,2);
M=N_1(:,1);
N=N_1(:,2);
%%% 调整C1和C2的坐标
[C1,L1]= newsecond00(C1,M)
[C2,L2]= newsecond00(C2,N)
C=[C1,C2]
%%% 根据maximin distance准则来调整坐标的大小
 [A , B]= newsecond01(C,N_1);
%%% 确定样本点的具体坐标值
A= newsecond02(A,L1,M);
B= newsecond02(B,L2,N);
C=[A,B]
N_1=[N_1;C];
    
end