function [ A ]= newsecond02( A ,L,N_1)
n=length(N_1);
m=length(A);
a=1/(2*(n+m));
%%%对新样本点进行插入排序与调整确定它的坐标
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
%%%%讨论两种初始点的插入情况
for k=1:m
    K=0;
 if d(1)<=a
     for i=2:n
         if l(i)>(2*(i-1)+1)*a
             t=l(i)-(2*(i-1)+1)*a;
             K=K+t;
             for j=i:n
                 l(j)=l(j)-t;end
             if A(k)<L
                 if A(k)<K
                     A(k)=(i-1)*2*a+A(k);
                     break;
                 end
             else
                 A(k)=N_1(n)+a+A(k)-L;
                 break;
             end
         end
     end
     if d(n+1)>a
         K=K+d(n+1)-a;
     end
 else
     K=d(1)-a;
     if  A(k)<K
         A(k)=A(k);
     else
          l(1)=a;
          for i=2:n
             l(i)=d(i)+l(i-1);end
          for i=2:n
               if l(i)>(2*(i-1)+1)*a
                  t=l(i)-(2*(i-1)+1)*a;
                  K=K+t;
                  for j=i:n
                      l(j)=l(j)-t;end
             if A(k)<L
                 if A(k)<K
                     A(k)=(i-1)*2*a+A(k);
                     break;
                 end
             else
                 A(k)=N_1(n)+a+A(k)-L;
                 break;
             end
             end
         end
     end      
 end
end
end