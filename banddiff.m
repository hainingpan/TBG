function banddiff(NV,Npsi,w)
Vrange=linspace(0,20,NV);
psirange=linspace(0,360,Npsi);
GammaMap=zeros(length(Vrange),length(psirange));
% kpMap=zeros(length(Vrange),length(psirange));
func=@(psi,V) energydiff(psi,V,w);

parfor Vindex=1:NV
    for psiindex=1:Npsi
        GammaMap(Vindex,psiindex)=feval(func,psirange(psiindex),Vrange(Vindex));
    end
end
save(sprintf('banddiffNV%dNpsi%dw%d.dat',NV,Npsi,w),'GammaMap','-ascii');

function [diffGamma,diffkp]=energydiff(psi,V,w)
parameters=mainTMD('psi',psi,'V',V,'w',w,'theta',3);
val=energyTMD(0,0,parameters);
diffGamma=val(1)-val(2);
% val=energyTMD(parameters.kp(1),parameters.kp(2),parameters);
% diffkp=val(1)-val(2);
% val=energyTMD(parameters.kn(1),parameters.kn(2),parameters);
% diffkn=val(1)-val(2);
end
end
