function [ N_1 ] = newfirst123( m,N_1 )
%NEWFIRST Summary of this function goes here
%   ��N_1�в���m����
n=length(N_1);
A=lhsamp(m,1);
%%% ��newsecond1������ȡA�ڷ�Ӱ�������е�����
[A,L]= newsecond1(A,N_1);
%%% ����Ӱ��������뵽������ƿռ���
A=newsecond2(A,L,N_1);
N_1 =[N_1;A];
[ N_1 ] = arrange( N_1 );
end