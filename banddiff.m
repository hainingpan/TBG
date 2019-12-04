Vrange=linspace(0,20,100);
psirange=linspace(0,360,100);
NVrange=length(Vrange);
Npsirange=length(psirange);
GammaMap=zeros(length(Vrange),length(psirange));
kpMap=zeros(length(Vrange),length(psirange));
knMap=zeros(length(Vrange),length(psirange));

parfor Vindex=1:NVrange
    for psiindex=1:Npsirange
        [GammaMap(Vindex,psiindex),kpMap(Vindex,psiindex),knMap(Vindex,psiindex)]=energydiff(psirange(psiindex),Vrange(Vindex),10);
    end
end
function [diffGamma,diffkp,diffkn]=energydiff(psi,V,w)
parameters=mainTMD('m',0.35,'psi',psi,'V',V,'w',w);
val=energyTMD(0,0,parameters);
diffGamma=val(1)-val(2);
val=energyTMD(parameters.kp(1),parameters.kp(2),parameters);
diffkp=val(1)-val(2);
val=energyTMD(parameters.kn(1),parameters.kn(2),parameters);
diffkn=val(1)-val(2);
end

