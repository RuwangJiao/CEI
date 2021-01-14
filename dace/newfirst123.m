function [ N_1 ] = newfirst123( m,N_1 )
%NEWFIRST Summary of this function goes here
%   在N_1中插入m个点
n=length(N_1);
A=lhsamp(m,1);
%%% 用newsecond1函数求取A在非影响区域中的坐标
[A,L]= newsecond1(A,N_1);
%%% 将非影响区域插入到整个设计空间中
A=newsecond2(A,L,N_1);
N_1 =[N_1;A];
[ N_1 ] = arrange( N_1 );
end