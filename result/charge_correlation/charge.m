a=zeros(12,1);
b=zeros(12,1);
for i=1:550
    if A(i,1)==1
        a(A(i,2)+8) = a(A(i,2)+8)+1;
    else
        b(A(i,2)+8) = b(A(i,2)+8)+1;
    end
end