[Neighborlist,tmat]=tb3_neighbor(neighborlist,t,3,parameters);

parameters.kpp=[-parameters.kp(1),parameters.kp(2)];
n=40;
counter=1;
N=3*n^2+3*n+1;
kxlist=zeros(N,1);
kylist=zeros(N,1);
kx3list=zeros(N,1);
ky3list=zeros(N,1);
a1=-(2*parameters.bM1+parameters.bM2)/3/n;
a2=(parameters.bM1+2*parameters.bM2)/3/n;
a31=a1/sqrt(3)*[0,1;-1,0]; %rotate 90 deg clockwisely
a32=a2/sqrt(3)*[0,1;-1,0];  %rotate 90 deg clockwisely
for yindex=-n:n
    for xindex=max(-n,-n+yindex):min(n+yindex,n)
        k=xindex*a1+yindex*a2;
        k2=xindex*a31+yindex*a32;
        kxlist(counter)=k(1);
        kylist(counter)=k(2);
        kx3list(counter)=k2(1);
        ky3list(counter)=k2(2);
        counter=counter+1;
    end
end
k=3;
% parfor i=1:N
%     entmp=energyTMD(kxlist(i)+parameters.bM2(1),kylist(i),parameters);
%     energyTMDlist(i)=max(entmp);
% end
energylist=real(tb([neighborlist{1:k+1}],[t{1:k+1}],kxlist,kylist,parameters));

energy3list=real(tb3(Neighborlist,tmat,kx3list,ky3list,parameters));


figure;
scatter(kxlist,kylist,[],energylist,'.');
