Vrange=linspace(0,20,50);
psirange=linspace(0,360,50);
bp=zeros(length(Vrange),length(psirange));
NV=length(Vrange);
Npsi=length(psirange);
parfor Vindex=1:NV
    for psiindex=1:Npsi
%         disp([Vindex,psiindex]);
        bp(Vindex,psiindex)=bf(psirange(psiindex),Vrange(Vindex),10);
    end
end        
save('bp.dat','bp','-ascii');


function bp=bf(psi,V,w)
parameters=mainTMD('m',0.35,'psi',psi,'V',V,'w',w,'theta',3);
[kcxmap,kcymap,kcx2map,kcy2map,bcmap,omega]=berrycurvature(1,parameters);
bp=sum(bcmap(:))*omega/(2*pi);
end
% save('berrycur.dat','bcmap','-ascii')
% save('kx.dat','kcx2map','-ascii')
% save('ky.dat','kcy2map','-ascii')