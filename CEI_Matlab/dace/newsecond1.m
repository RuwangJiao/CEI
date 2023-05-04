function [A,L] = newsecond1(A,N_1)
%NEWSECOND Summary of this function goes here
% 确定x,y上的坐标 
n=length(N_1);
m=length(A);
a=1/(2*(n+m));
%%%对新样本点进行插入排序与调整确定它的坐标
[ N_1 ] = arrange( N_1 );
d(1)=N_1(1);
d(n+1)=1-N_1(n);
l(1)=d(1);
K=0;
for k=2:n
    d(k)=N_1(k)-N_1(k-1);
end
for i=2:n
    l(i)=d(i)+l(i-1);
end
%%%确定K的大小来求解A的值
 if d(1)<=a
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
 K
 A=A*K
 %%%  对A中的点据maximin distance 准则进行调整
A=arrange(A);
for i=1:m-1
    if A(i+1)-A(i)<1/(2*(n+m))
        A(i+1)=A(i)+1/(2*(n+m));
    end
end


