function [ C1,C2 ]= newsecond01( D ,N_1)
 %%%  对A中的点据maximin distance 准则进行调整N_1是初始样本点D是新增的点
 A1=D(:,1);
 A2=D(:,2);
 m=length(A1);
 n=length(N_1(:,1));
 B1=arrange(A1)
 for i=1:m
     for j=1:m
         if B1(i,1)==A1(j)
             B2(i,1)=A2(j)
             break;
         end
     end
 end
for i=1:m-1
    if 1/(2*(n+m))<B1(i+1)- B1(i)
        B1(i+1)=B1(i)+1/(2*(n+m));
    end
end
 C2=arrange(B2);
 for i=1:m
     for j=1:m
         if C2(i)==B2(j)
             C1(i,1)=B1(j);
         end
     end
 end
for i=1:m-1
    if 1/(2*(n+m))<C2(i+1)-C2(i)
        C2(i+1)=C2(i)+1/(2*(n+m));
    end
end
end