function [A,L] = newsecond00(A,N_1)
%NEWSECOND Summary of this function goes here
% ȷ��x,y�ϵ�����
n=length(N_1);
m=length(A);
a=1/(2*(n+m));
%%%���A������
[ N_1 ] = arrange( N_1 );
d(1)=N_1(1);
d(n+1)=1-N_1(n);
l(1)=d(1);
for k=2:n
    d(k)=N_1(k)-N_1(k-1);
end
for i=2:n
    l(i)=d(i)+l(i-1);
end
%%%ȷ��K�Ĵ�С�����A��ֵ
 if d(1)<=a
     K=0;
     for i=2:n
         if l(i)>(2*(i-1)+1)*a
             t=l(i)-(2*(i-1)+1)*a;
             K=K+t;
             for j=i:n
                 l(j)=l(j)-t;end
         end
     end
          L=K;
     if d(n+1)>a
         K=K+d(n+1)-a;
     end
 else
          K=d(1)-a;
          l(1)=a;
          for i=2:n
             l(i)=d(i)+l(i-1);end
          for i=2:n
               if l(i)>(2*(i-1)+1)*a
                  t=l(i)-(2*(i-1)+1)*a;
                  K=K+t;
                  for j=i:n
                      l(j)=l(j)-t;end
             end
          end 
          L=K;
      if d(n+1)>a
         K=K+d(n+1)-a;
     end
 end
 A=A*K;
 K
end