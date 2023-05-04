function [R_square RAAE RMAE]=metrics(y1,y2)
% ���������ں��������ܲ��ԣ��������������ݵ����
% y1����ʵ���ݣ�y2��Ԥ������
m=size(y1,1);
STD=std(y1);

% ��������ָ��R_square
t=0;
k=0;
for i_test=1:m
    t=t+(y1(i_test)-y2(i_test))^2;
     k=k+(y1(i_test)-mean(y2))^2;
end
R_square=1-t/k;

% ��������ָ��RAAE
t=0;
for i=1:m
    t=t+abs(y1(i)-y2(i));
end
RAAE=t/(m*STD);

% ��������ָ��RMAE
k=zeros(1,m);
for i=1:m
    k(i)=abs(y1(i)-y2(i));
end
RMAE=max(k)/STD;

end