function run_dos_vs_nu()
Vzlist=-28:1:-28;
n=20;
N=3*n^2+3*n+1;
Nen=2000;
kxlist=zeros(1,N);
kylist=zeros(1,N);
energylist=zeros(length(Vzlist),N);
enlist=zeros(length(Vzlist),Nen);
nu=zeros(length(Vzlist),Nen);
dos=zeros(length(Vzlist),Nen);
E_vanHove=zeros(1,length(Vzlist));
gap=zeros(1,length(Vzlist));
for i=1:length(Vzlist)
    if i==1
        [kxlist(i,:),kylist(i,:),energylist(i,:),enlist(i,:),nu(i,:),dos(i,:),E_vanHove(i)]=dos_vs_nu(4,Vzlist(i),700,Nen);
    else
        [~,~,energylist(i,:),enlist(i,:),nu(i,:),dos(i,:),E_vanHove(i)]=dos_vs_nu(4,Vzlist(i),700,Nen);
    end
    [gap(i),~]=phasediagram(4,Vzlist(i),0.02);
end
clear i n N Nen
save('dosvsnuVz-28.mat');
% save('gap.dat','gap','-ascii');
end