function LDOS_TMD_phase(NV,Npsi,w)
Vrange=linspace(0,20,NV);
psirange=linspace(0,360,Npsi);
dos1=zeros(length(Vrange),length(psirange),50);
dos2=zeros(length(Vrange),length(psirange),50);
enlist1=zeros(length(Vrange),length(psirange),50);
enlist2=zeros(length(Vrange),length(psirange),50);
% NV=length(Vrange);
% Npsi=length(psirange);
parfor Vindex=1:NV
    for psiindex=1:Npsi
%         disp([Vindex,psiindex]);
        psi=psirange(psiindex);
        V=Vrange(Vindex);
        parameters=mainTMD('m',0.35,'psi',psi,'V',V,'w',w,'theta',3);
        [dos1(Vindex,psiindex,:),dos2(Vindex,psiindex,:),...
            enlist1(Vindex,psiindex,:),enlist2(Vindex,psiindex,:)]=LDOS_TMD_r(parameters);
    end
end

save(sprintf('dosNV%dNpsi%dw%d.mat',NV,Npsi,w));
end
% parameters=mainTMD('m',0.35,'psi',-0.3329/(2*pi)*360,'V',4.428,'w',20,'theta',3);
% parameters=mainTMD('m',0.35,'psi',180,'V',0,'w',10,'theta',3);

% parameters=mainTMD('m',0.35,'psi',100,'V',20,'w',10);
% [dos1,dos2,enlist1,enlist2]=LDOS_TMD_r(parameters);
% [int1,int2]=LDOS_TMD_r(parameters);

% [ldos,enlist]=LDOS_TMD_rx(parameters);