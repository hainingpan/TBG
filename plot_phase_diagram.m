Nangle=50;
NVz=50;
anglelist=linspace(3.5,5.5,Nangle);
Vzlist=linspace(0,50,NVz);
epsilonlist=0.005:0.005:0.04;
gap=zeros(NVz,Nangle,length(epsilonlist));
isconverge=zeros(NVz,Nangle,length(epsilonlist));

for i=1:NVz
    for j=1:Nangle
        fprintf()
        [gap(i,j,:),isconverge(i,j,:)]=phasediagram(anglelist(j),Vzlist(i),epsilonlist);
    end
end

save('phasediagram.mat','gap','isconverge','anglelist','Vzlist','epsilonlist');