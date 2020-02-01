NVz=50;
Vzlist=linspace(0,50,NVz);
Neps=50;
epsilonlist=linspace(0.005,0.04,Neps);
thetalist=[4,4.2,5.1];
gap=zeros(NVz,length(epsilonlist),length(thetalist));
isconverge=zeros(NVz,length(epsilonlist),length(thetalist));

for k=1:3
    parfor i=1:NVz
            [gap(i,:,k),isconverge(i,:,k)]=phasediagram(thetalist(k),Vzlist(i),epsilonlist);
    end
end
