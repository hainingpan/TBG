Vrange=linspace(0,20,100);
psirange=linspace(0,360,100);
bp=zeros(length(Vrange),length(psirange));

for Vindex=1:length(Vrange)
    for psiindex=1:length(psirange)
        disp([Vindex,psiindex]);
        bp(Vindex,psiindex)=bf(psirange(psiindex),Vrange(Vindex),10);
    end
end

        



function bp=bf(psi,V,w)
parameters=mainTMD('m',0.35,'psi',psi,'V',V,'w',w);
[kcxmap,kcymap,kcx2map,kcy2map,bcmap,omega]=berrycurvature(1,parameters);
bp=sum(bcmap(:))*omega/(2*pi);
end
% save('berrycur.dat','bcmap','-ascii')
% save('kx.dat','kcx2map','-ascii')
% save('ky.dat','kcy2map','-ascii')