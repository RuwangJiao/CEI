function S=inheritlhd(m,n)
%LHD Summary of this function goes here
S=lhsamp(m,2);
a=0.005;
for i=1:n
    A=lhsamp(m,2);
    for k=1:m-1
         A(k,:)=A(k,:)*(1-2*a*m);
        if A(k,1)<S(k,1)
             A(k,1)=A(k,1);
        elseif S(k,1)<A(k,1)<S(k+1,1)
            A(k,1)=A(k,1)+2*k*a;
        else  A(k,1)=A(k,1)+2*m*a;end
        if A(k,2)<S(k,2)
             A(k,2)=A(k,2);
        elseif S(k,2)<A(k,2)<S(k+1,2)
            A(k,2)=A(k,2)+2*k*a;
        else  A(k,2)=A(k,2)+2*m*a;end
     end
A
 S=[S;A];
end
S

