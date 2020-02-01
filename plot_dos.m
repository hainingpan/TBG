Nangle=50;
NVz=100;
anglelist=linspace(3.,6,Nangle);
Vzlist=linspace(0,50,NVz);
energy_at_half=zeros(NVz,Nangle);
dos=zeros(NVz,Nangle);

parfor i=1:NVz
    for j=1:Nangle
        [energy_at_half(i,j),dos(i,j)]=dos_at_half_filling(anglelist(j),Vzlist(i));
    end
end

save('dos.mat','anglelist','Vzlist','energy_at_half','dos');