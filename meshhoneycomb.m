n=20;
a1=[sqrt(3)/2,-1/2]/n;
a2=[0,1]/n;
counter=1;
k=zeros(3*n^2+3*n+1,2);
for yindex=-n:n
    for xindex=max(-n,-n+yindex):min(n+yindex,n)
        k(counter,:)=xindex*a1+yindex*a2;
        counter=counter+1;
    end
end
        
        