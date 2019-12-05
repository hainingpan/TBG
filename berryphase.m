function berryphase(NV,Npsi,w)
Vrange=linspace(0,20,50);
psirange=linspace(0,360,50);
bp=zeros(length(Vrange),length(psirange));
% NV=length(Vrange);
% Npsi=length(psirange);
func=@(psi,V) bf(psi,V,w);
parfor Vindex=1:NV
    for psiindex=1:Npsi
%         disp([Vindex,psiindex]);
        bp(Vindex,psiindex)=feval(func,psirange(psiindex),Vrange(Vindex));
    end
end        
save(sprintf('bpNV%dNpsi%dw%d.dat',NV,Npsi,w),'bp','-ascii');

function bp=bf(psi,V,w)
parameters=mainTMD('m',0.35,'psi',psi,'V',V,'w',w,'theta',3);
[~,~,~,~,bcmap,omega]=berrycurvature(1,parameters);
bp=sum(bcmap(:))*omega/(2*pi);
end
end