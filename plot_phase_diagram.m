Nangle=50;
NVz=50;
anglelist=linspace(3.,6,Nangle);
% anglelist=[3,4];
% Nangle=length(anglelist);
Vzlist=linspace(0,50,NVz);
% epsilonlist=0.005:0.005:0.04;
% epsilonlist=1./[10,20,30,40];
epsilonlist=1./[30];
gap=zeros(NVz,Nangle,length(epsilonlist));
isconverge=zeros(NVz,Nangle,length(epsilonlist));

parfor i=1:NVz
    for j=1:Nangle
        [gap(i,j,:),isconverge(i,j,:)]=phasediagram(anglelist(j),Vzlist(i),epsilonlist);
    end
end

save('phasediagram.mat','gap','isconverge','anglelist','Vzlist','epsilonlist');