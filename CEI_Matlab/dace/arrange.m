function [ R_1 ] = arrange( N_1 )
%ARRANGE Summary of this function goes here
%  对现有的采样点从小到大排序
i=size(N_1,1); 
for j=1:i
     for k=j:i
         if N_1(j)>N_1(k)
             t=N_1(j);
             N_1(j)=N_1(k);
             N_1(k)=t;
         else N_1(j)=N_1(j);
         end
     end
end
 R_1=N_1;


end

