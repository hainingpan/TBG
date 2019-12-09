function re=hoppingt(an,parameters)
n=5;
state=1;
xrange=-n:n;
yrange=-n:n;
bM1=parameters.bM1;
bM2=parameters.bM2;


%shift to diamond
% a1=-bM1/(2*n);
% a2=(bM1+bM2)/(2*n);
%shift to rectangular
a1=bM2/(2*n);
a2=(-2*bM1-bM2)/2/(2*n);
enmap=zeros(2*n+1,2*n+1);
kxmap=zeros(2*n+1,2*n+1);
kymap=zeros(2*n+1,2*n+1);


Nx=length(xrange);
Ny=length(yrange);
omega=abs(cross([bM1,0],[bM2,0]));
omega=omega(3);
parfor xindex=1:Nx
    kx=xrange(xindex);
    for yindex=1:Ny
        ky=yrange(yindex);
        k=kx*a1+ky*a2;
        [val,~]=energyTMD(k(1),k(2),parameters);
        enmap(xindex,yindex)=val(state);
        kxmap(xindex,yindex)=k(1);
        kymap(xindex,yindex)=k(2);
    end
end

for i=1:length(an)
    neigh=an{i}(1)*parameters.aM1+an{i}(2)*parameters.aM2;

    %use summation
    %re(i)=sum(exp(1i*(kxmap*neigh(1)+kymap*neigh(2))).*enmap,'all')/(2*n+1)^2;

    %use trapzodial integral
    F=exp(1i*(kxmap*neigh(1)+kymap*neigh(2))).*enmap;
    re(i)=trapz(kymap(1,:),trapz(kxmap(:,1),F,2))/omega;
end
end
